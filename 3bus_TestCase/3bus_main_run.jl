#!/usr/bin/env julia

#Main script for balanced EMT transient simulation of 3 with SG, GFM, and GFL resources.
#also Includes additinal script for  small-signal stability analysis.

using PowerSystems
using PowerSimulationsDynamics
using PowerFlows
using Sundials
using OrdinaryDiffEq
using Plots
using LinearAlgebra
using Statistics
using Dates
using Printf
using Logging
using CSV
using DataFrames

const PSY  = PowerSystems
const PSID = PowerSimulationsDynamics

include("case3_matpower.jl")
include("3bus_parameters.jl")

gr()
default(fontfamily = "sans-serif")

#Run event: load step
const LOAD_CHANGE_ENABLE = true      #  Enable load-step event
const LOAD_STEP_BUS      = 2   #Load-step bus bus2
const LOAD_STEP_P_MW     = 25.0     #Active-power step[MW]
const LOAD_STEP_Q_MVAr   = 5.0     #Reactive-power step[MVAr]
const LOAD_STEP_T0       = 0.25   #Load-step time[s]

const LOAD_CHANGE_EVENTS = [
    Dict(
        :buses        => [LOAD_STEP_BUS],#Buses with load step
        :event_time_s => LOAD_STEP_T0,   #Event time[s]
        :p_pu_map     => Dict(LOAD_STEP_BUS => LOAD_STEP_P_MW / SYS_BASE_MVA),   #P target[p.u.]
        :q_pu_map     => Dict(LOAD_STEP_BUS => LOAD_STEP_Q_MVAr / SYS_BASE_MVA),     #Q target[p.u.]
    ),
]

#Labels
function model_label(b::Int)
    if b in SG_BUSES
        return b == SLACK_BUS_NUM ? "SG-slack" : "SG"
    elseif b in GFM_BUSES
        return b == SLACK_BUS_NUM ? "GFM-slack" : "GFM"
    elseif b in GFL_BUSES
        return "GFL"
    end
    return "UNSET"
end

bus_label(b::Int) = @sprintf("Bus %d (%s)", b, model_label(b))

function add_component_clean!(args...)
    with_logger(NullLogger()) do
        add_component!(args...)
    end
end

#helper functions
safe_pinj(pf, idx) = try pf.bus_activepower_injection[idx]   catch; 0.0 end
safe_qinj(pf, idx) = try pf.bus_reactivepower_injection[idx] catch; 0.0 end
safe_pwd(pf, idx)  = try pf.bus_activepower_withdrawals[idx] catch; 0.0 end
safe_qwd(pf, idx)  = try pf.bus_reactivepower_withdrawals[idx] catch; 0.0 end
safe_gen_base(g::Generator) = try PSY.get_base_power(g) catch; SYS_BASE_MVA end

#Validation
function validate_partitions!(sg::Vector{Int}, gfm::Vector{Int}, gfl::Vector{Int})
    sgset, gfmset, gflset = Set(sg), Set(gfm), Set(gfl)
    allowed = Set(1:3)

    issubset(sgset, allowed)  || error("SG_BUSES must be within 1:3.")
    issubset(gfmset, allowed) || error("GFM_BUSES must be within 1:3.")
    issubset(gflset, allowed) || error("GFL_BUSES must be within 1:3.")

    ov1, ov2, ov3 = intersect(sgset, gfmset), intersect(sgset, gflset), intersect(gfmset, gflset)
    (isempty(ov1) && isempty(ov2) && isempty(ov3)) || error("Partitions overlap.")

    union(sgset, gfmset, gflset) == allowed || error("Partitions must cover buses 1:3 exactly.")
    SLACK_BUS_NUM in gflset && error("Slack bus cannot be GFL; use SG or GFM for the reference bus.")
    return true
end

#Load-event schedule
struct RealizedLoadStep
    bus::Int
    t_step::Float64
    ΔP_MW::Float64
    ΔQ_MVAr::Float64
    refP_MW::Float64
    refQ_MVAr::Float64
    refP_pu::Float64
    refQ_pu::Float64
    deltaP_pu::Float64
    deltaQ_pu::Float64
end

_event_p_map(ev) = haskey(ev, :p_pu_map) ? ev[:p_pu_map] : get(ev, :pu_map, Dict{Int,Float64}())
_event_q_map(ev) = get(ev, :q_pu_map, Dict{Int,Float64}())
_to_int_float(raw) = Dict{Int,Float64}(Int(k) => Float64(v) for (k, v) in raw)

function build_realized_load_steps()
    prevP = Dict{Int,Float64}()
    prevQ = Dict{Int,Float64}()
    steps = RealizedLoadStep[]
    for ev in sort(LOAD_CHANGE_EVENTS, by = e -> Float64(e[:event_time_s]))
        t = Float64(ev[:event_time_s])
        p_map = _to_int_float(_event_p_map(ev))
        q_map = _to_int_float(_event_q_map(ev))
        for b in Vector{Int}(ev[:buses])
            refP_pu = get(p_map, b, 0.0); refQ_pu = get(q_map, b, 0.0)
            refP_MW = refP_pu * SYS_BASE_MVA; refQ_MVAr = refQ_pu * SYS_BASE_MVA
            ΔP = refP_MW - get(prevP, b, 0.0); ΔQ = refQ_MVAr - get(prevQ, b, 0.0)
            prevP[b] = refP_MW; prevQ[b] = refQ_MVAr
            push!(steps, RealizedLoadStep(b, t, ΔP, ΔQ, refP_MW, refQ_MVAr,
                                          refP_pu, refQ_pu, ΔP / SYS_BASE_MVA, ΔQ / SYS_BASE_MVA))
        end
    end
    sort!(steps, by = s -> s.t_step)
    return steps
end

function print_load_schedule(steps::Vector{RealizedLoadStep})
    println("\n" * "="^96)
    println("ACTIVE/REACTIVE LOAD-STEP EVENT SCHEDULE")
    println("="^96)
    @printf("%-8s %-10s %-22s %-22s %-22s\n", "Bus", "t [s]", "ΔP*", "ΔQ*", "P*/Q* ref")
    println("-"^96)
    for s in steps
        ps = s.ΔP_MW >= 0 ? "+" : ""; qs = s.ΔQ_MVAr >= 0 ? "+" : ""
        @printf("%-8d %-10.3f %-22s %-22s %-22s\n", s.bus, s.t_step,
                @sprintf("%s%.1f MW", ps, s.ΔP_MW), @sprintf("%s%.1f MVAr", qs, s.ΔQ_MVAr),
                @sprintf("%.1f MW / %.1f MVAr", s.refP_MW, s.refQ_MVAr))
    end
    println("="^96 * "\n")
end

function event_summary(steps::Vector{RealizedLoadStep})
    isempty(steps) && return "No load events"
    return join([@sprintf("Bus %d, t=%.2fs: ΔP=%+.1f MW, ΔQ=%+.1f MVAr",
                          s.bus, s.t_step, s.ΔP_MW, s.ΔQ_MVAr) for s in steps], "  |  ")
end

#System utilities
function gen_at_bus(sys::System, busnum::Int)
    for g in get_components(Generator, sys)
        PSY.get_number(PSY.get_bus(g)) == busnum && return g
    end
    error("No Generator found at bus $busnum")
end

function bus_component(sys::System, busnum::Int)
    for b in get_components(Bus, sys)
        PSY.get_number(b) == busnum && return b
    end
    error("Bus $busnum not found")
end

function update_system_voltages!(sys::System, pf::PowerFlowData)
    for bus in get_components(Bus, sys)
        bn = PSY.get_number(bus)
        haskey(pf.bus_lookup, bn) || continue
        idx = pf.bus_lookup[bn]
        PSY.set_magnitude!(bus, pf.bus_magnitude[idx])
        PSY.set_angle!(bus, pf.bus_angles[idx])
    end
end

#Use\PF generator injections so device setpoints match the solved case.
function sync_static_gens_to_pf!(sys::System, pf::PowerFlowData; buses::Vector{Int} = collect(1:3))
    for b in buses
        haskey(pf.bus_lookup, b) || continue
        idx = pf.bus_lookup[b]
        g = gen_at_bus(sys, b)
        PSY.set_active_power!(g, safe_pinj(pf, idx))
        try PSY.set_reactive_power!(g, safe_qinj(pf, idx)) catch end
    end
end

function enforce_load_bases!(sys::System)
    for ld in get_components(ElectricLoad, sys)
        try PSY.set_base_power!(ld, SYS_BASE_MVA) catch end
    end
end

function set_generator_bases!(sys::System)
    for g in get_components(Generator, sys)
        b = PSY.get_number(PSY.get_bus(g))
        try PSY.set_base_power!(g, getGenMBase(b)) catch e
            @warn "Could not set generator base power" bus = b err = string(e)
        end
    end
end

#Load model conversion
function loads_to_constant_impedance!(sys::System)
    converted = 0
    for ld in collect(get_components(PowerLoad, sys))
        nm = PSY.get_name(ld)
        startswith(nm, "PerturbLoad") && continue
        b = PSY.get_bus(ld)
        V = max(try PSY.get_magnitude(b) catch; 1.0 end, 0.5)
        P = try PSY.get_active_power(ld)   catch; 0.0 end
        Q = try PSY.get_reactive_power(ld) catch; 0.0 end
        (abs(P) < 1e-9 && abs(Q) < 1e-9) && continue
        Y = complex(P, -Q) / V^2
        adm = nothing
        try
            adm = FixedAdmittance(name = "Zload_" * nm, available = true, bus = b, Y = Y)
        catch
            adm = FixedAdmittance("Zload_" * nm, true, b, Y)
        end
        remove_component!(sys, ld)
        add_component!(sys, adm)
        converted += 1
    end
    println("Converted $converted constant-power loads to constant impedance (FixedAdmittance).")
    return converted
end

# Dynamic-line helpers
normalize_pair(i::Int, j::Int) = i < j ? (i, j) : (j, i)

function line_pair(line)
    try
        arc = PSY.get_arc(line)
        return normalize_pair(Int(PSY.get_number(arc.from)), Int(PSY.get_number(arc.to)))
    catch
        buses = PSY.get_buses(line)
        return normalize_pair(Int(PSY.get_number(buses[1])), Int(PSY.get_number(buses[2])))
    end
end

available_line_pairs(sys::System) = unique(sort([line_pair(l) for l in get_components(Line, sys)]))

function selected_dynamic_line_pair_set(sys::System)
    allpairs = Set{Tuple{Int,Int}}(available_line_pairs(sys))
    forced_static = Set(normalize_pair(p[1], p[2]) for p in DYNAMIC_LINE_FORCE_STATIC)
    USE_SELECTIVE_DYNAMIC_LINES || return setdiff(allpairs, forced_static)
    chosen = Set{Tuple{Int,Int}}()
    for p in DYNAMIC_LINE_INCLUDE
        np = normalize_pair(p[1], p[2])
        np in allpairs ? push!(chosen, np) :
            @warn "Requested dynamic line not found" pair = np available = sort(collect(allpairs))
    end
    return setdiff(chosen, forced_static)
end

function print_dynamic_line_selection(sys::System, include_pairs::Set{Tuple{Int,Int}})
    allpairs = Set{Tuple{Int,Int}}(available_line_pairs(sys))
    println("\n" * "="^90)
    println("DYNAMIC-LINE CONFIGURATION")
    println("="^90)
    println("Available AC lines = $(sort(collect(allpairs)))")
    println("Dynamic lines      = $(sort(collect(include_pairs)))")
    println("Static lines       = $(sort(collect(setdiff(allpairs, include_pairs))))")
    println("="^90 * "\n")
end

function apply_dynamic_line_damping_edits!(sys::System, include_pairs::Set{Tuple{Int,Int}})
    rows = NamedTuple[]
    APPLY_DYNAMIC_LINE_DAMPING_EDITS || return rows
    for line in get_components(Line, sys)
        pair = line_pair(line)
        pair in include_pairs || continue
        r_old = NaN; x_old = NaN
        try
            r_old = PSY.get_r(line); x_old = PSY.get_x(line)
        catch
            @warn "Could not read line R/X for damping edits" line = PSY.get_name(line) pair = pair
            continue
        end
        (isfinite(r_old) && isfinite(x_old)) || continue
        r_new = r_old <= 0.0 ? DYNAMIC_LINE_R_FLOOR_PU : r_old * DYNAMIC_LINE_R_SCALE
        x_new = x_old * DYNAMIC_LINE_X_SCALE
        try
            PSY.set_r!(line, r_new); PSY.set_x!(line, x_new)
        catch err
            @warn "Could not apply R/X edit" line = PSY.get_name(line) pair = pair err = string(err)
            continue
        end
        push!(rows, (line = PSY.get_name(line), pair = pair, r_old = r_old, r_new = r_new, x_old = x_old, x_new = x_new))
    end
    return rows
end

function print_dynamic_line_damping_edits(rows)
    isempty(rows) && return
    println("\n" * "="^96)
    println("DYNAMIC-LINE SERIES R/X EDITS")
    println("="^96)
    @printf("%-22s %-12s %-12s %-12s %-12s %-12s\n", "Line", "Pair", "R_old", "R_new", "X_old", "X_new")
    println("-"^96)
    for r in rows
        @printf("%-22s %-12s %-12.6f %-12.6f %-12.6f %-12.6f\n",
                r.line, string(r.pair), r.r_old, r.r_new, r.x_old, r.x_new)
    end
    println("="^96 * "\n")
end

#Add selected dynamic branches.
function add_dynamic_lines!(sys::System; include_pairs::Set{Tuple{Int,Int}} = selected_dynamic_line_pair_set(sys))
    added = Tuple{Int,Int}[]; skipped = Tuple{Int,Int}[]
    for line in collect(get_components(Line, sys))
        pair = line_pair(line)
        if !(pair in include_pairs)
            push!(skipped, pair); continue
        end
        try
            add_component_clean!(sys, DynamicBranch(line)); push!(added, pair)
        catch e
            @warn "DynamicBranch skipped" line = PSY.get_name(line) pair = pair err = string(e)
        end
    end
    println("DynamicBranch added = $(length(unique(added))): $(sort(unique(added)))")
    return unique(sort(added)), unique(sort(skipped))
end

# Reference helpers
function get_pref_qref(g::Generator)
    p = try PSY.get_active_power(g)   catch; 0.0 end
    q = try PSY.get_reactive_power(g) catch; 0.0 end
    scale = SYS_BASE_MVA / safe_gen_base(g)
    return p * scale, q * scale
end

get_vref(sys::System, g::Generator) = try PSY.get_magnitude(PSY.get_bus(g)) catch; 1.0 end

function print_reference_summary(sys::System)
    println("\n" * "="^104)
    println("DYNAMIC MODEL REFERENCE SUMMARY (device base, from PF operating point)")
    println("="^104)
    @printf("%-8s %-10s %-16s %-16s %-14s %-12s %-12s\n",
            "Bus", "Model", "P_ref[pu,dev]", "Q_ref[pu,dev]", "V_ref[pu]", "mBase[MVA]", "P_ref[MW]")
    println("-"^104)
    for b in 1:3
        g = gen_at_bus(sys, b)
        p, q = get_pref_qref(g)
        mb = safe_gen_base(g)
        @printf("%-8d %-10s %-16.6f %-16.6f %-14.6f %-12.2f %-12.3f\n",
                b, model_label(b), p, q, get_vref(sys, g), mb, p * mb)
    end
    println("="^104 * "\n")
end

#Perturbable loads
function ensure_perturb_load!(sys::System, busnum::Int)
    name = @sprintf("PerturbLoad_bus%02d", busnum)
    for ld in get_components(ElectricLoad, sys)
        if PSY.get_name(ld) == name
            try PSY.set_base_power!(ld, SYS_BASE_MVA) catch end
            return ld
        end
    end
    b = bus_component(sys, busnum)
    ld = nothing
    try
        ld = PowerLoad(name = name, available = true, bus = b, active_power = 0.0,
                       reactive_power = 0.0, base_power = SYS_BASE_MVA,
                       max_active_power = 20.0, max_reactive_power = 20.0)
    catch
        try
            ld = PowerLoad(name, true, b, 0.0, 0.0, SYS_BASE_MVA, 20.0, 20.0)
        catch
            ld = PowerLoad(name, true, b, 0.0, 0.0, SYS_BASE_MVA)
        end
    end
    add_component!(sys, ld)
    try PSY.set_base_power!(ld, SYS_BASE_MVA) catch end
    return ld
end

function perturb_load_at_bus(sys::System, busnum::Int)
    target = @sprintf("PerturbLoad_bus%02d", busnum)
    for ld in get_components(ElectricLoad, sys)
        PSY.get_name(ld) == target && return ld
    end
    error("Perturb load $target not found")
end

function make_load_change(sys::System; bus::Int, t_step::Float64, setval::Float64, channel::Symbol)
    ld = perturb_load_at_bus(sys, bus)
    set_pu = setval / SYS_BASE_MVA
    fields = channel === :P ? (:active_power, :P, :p, :P_ref_power, :P_ref) :
                              (:reactive_power, :Q, :q, :Q_ref_power, :Q_ref)
    for field in fields
        try return LoadChange(t_step, ld, field, set_pu) catch end
    end
    error("Could not apply $(channel) load step on $(typeof(ld)).")
end

function make_load_changes_from_event_schedule(sys::System, steps::Vector{RealizedLoadStep})
    perts = PSID.Perturbation[]
    for s in steps

# LoadChange uses absolute P/Q targets.
        abs(s.ΔP_MW)   > 1e-12 && push!(perts, make_load_change(sys; bus = s.bus, t_step = s.t_step, setval = s.refP_MW,   channel = :P))
        abs(s.ΔQ_MVAr) > 1e-12 && push!(perts, make_load_change(sys; bus = s.bus, t_step = s.t_step, setval = s.refQ_MVAr, channel = :Q))
    end
    return perts
end

#Dynamic model attachment
sg_machine_params() = Dict{Symbol,Float64}(
    :Rs => 0.02, :Xl => 0.300, :Xd => 2.000, :Xq => 1.900,
    :Xd_p => 0.600, :Xq_p => 0.800, :Xd_pp => 0.400,
    :Tdo_p => 4.200, :Tqo_p => 0.700, :Tdo_pp => 0.040, :Tqo_pp => 0.035)

sg_avr_params() = Dict{Symbol,Float64}(
    :Ka => 20.0, :Ke => 1.0, :Kf => 0.1, :Ta => 0.02, :Te => 0.5, :Tf => 1.0, :Tr => 0.05)

function make_vi_block(; Ta::Float64, kd::Float64, kω::Float64, P_ref::Float64)
    try return VirtualInertia(Ta = Ta, kd = kd, kω = kω, P_ref = P_ref)
    catch; return VirtualInertia(Ta = Ta, kd = kd, k_w = kω, P_ref = P_ref) end
end

function attach_sg!(sys::System, g::Generator;
        H::Float64, D::Float64 = D_SG_DEFAULT, droop::Float64 = R_SG_PU,
        mach_p::Dict{Symbol,Float64}, avr_p::Dict{Symbol,Float64})

    machine = RoundRotorQuadratic(R = mach_p[:Rs], Td0_p = mach_p[:Tdo_p], Td0_pp = mach_p[:Tdo_pp],
        Tq0_p = mach_p[:Tqo_p], Tq0_pp = mach_p[:Tqo_pp], Xd = mach_p[:Xd], Xq = mach_p[:Xq],
        Xd_p = mach_p[:Xd_p], Xq_p = mach_p[:Xq_p], Xd_pp = mach_p[:Xd_pp], Xl = mach_p[:Xl], Se = (0.0, 0.0))
    shaft = SingleMass(H = H, D = D)
    avr = AVRTypeI(Ka = avr_p[:Ka], Ke = avr_p[:Ke], Kf = avr_p[:Kf], Ta = avr_p[:Ta],
        Te = avr_p[:Te], Tf = avr_p[:Tf], Tr = avr_p[:Tr], Va_lim = (-10.0, 10.0), Ae = 0.001, Be = 1.0)
    tg = TGTypeI(R = droop, Ts = TG_TS, Tc = TG_EPS, T3 = TG_EPS, T4 = TG_EPS, T5 = TG_EPS,
        valve_position_limits = (min = -5.0, max = 5.0))

# Use stateless PSS when stabilizer gain is zero.
    pss = PSSFixed(V_pss = 0.0)
    dyn = DynamicGenerator(name = PSY.get_name(g), ω_ref = 1.0, machine = machine,
        shaft = shaft, avr = avr, prime_mover = tg, pss = pss)
    add_component_clean!(sys, dyn, g)
    return dyn
end

# GFM filter and inner-loop 
function attach_gfm!(sys::System, g::Generator;
        droop::Float64 = R_GFM_PU, kq::Float64 = 0.03, Ta::Float64 = GFM_TA_DEFAULT,
        kd::Float64 = GFM_KD_DEFAULT, kω::Float64 = GFM_KW_DEFAULT,
        kpv::Float64 = 0.8, kiv::Float64 = 300.0, kpc::Float64 = 1.0, kic::Float64 = 120.0,
        rv::Float64 = 0.0, lv::Float64 = 0.0, ωad::Float64 = 1500.0, kad::Float64 = 20.0,
        lf::Float64 = 0.08, rf::Float64 = 0.030, cf::Float64 = 0.074, lg::Float64 = 0.2, rg::Float64 = 0.05)

    p_ref0, q_ref0 = get_pref_qref(g)
    v_ref0 = get_vref(sys, g)
    converter = AverageConverter(rated_voltage = 1.0, rated_current = 1.0)
    vi = make_vi_block(Ta = Ta, kd = kd, kω = kω, P_ref = p_ref0)
    outer = OuterControl(vi, ReactivePowerDroop(kq = kq, ωf = 10.0, V_ref = v_ref0))
    inner = VoltageModeControl(kpv = kpv, kiv = kiv, kffv = 0.0, rv = rv, lv = lv,
        kpc = kpc, kic = kic, kffi = 0.0, ωad = ωad, kad = kad)
    dc  = FixedDCSource(voltage = 1200.0)
    pll = FixedFrequency()
    filt = LCLFilter(lf = lf, rf = rf, cf = cf, lg = lg, rg = rg)
    dyn = DynamicInverter(name = PSY.get_name(g), ω_ref = 1.0, converter = converter,
        outer_control = outer, inner_control = inner, dc_source = dc, freq_estimator = pll, filter = filt)
    add_component_clean!(sys, dyn, g)
    return dyn
end

#GFL outer loop and LCL 
function attach_gfl!(sys::System, g::Generator;
        kp_p::Float64 = 0.35,
        ki_p::Float64 = 4.0,
        ωz_p::Float64 = 40.0,
        kp_q::Float64 = 0.35,
        ki_q::Float64 = 4.0,
        ωf_q::Float64 = 40.0,
        kω_droop::Float64 = GFL_KW_DEFAULT,
        kpc::Float64 = 0.7,
        kic::Float64 = 60.0,
        kffv::Float64 = 0.0,
        pll_ωlp::Float64 = 50.0,
        pll_kp::Float64 = GFL_PLL_KP,
        pll_ki::Float64 = GFL_PLL_KI,
        lf::Float64 = 0.08,
        rf::Float64 = 0.04,
        cf::Float64 = 0.074,
        lg::Float64 = 0.2,
        rg::Float64 = 0.06)

    p_ref0, q_ref0 = get_pref_qref(g)
    converter = AverageConverter(rated_voltage = 1.0, rated_current = 1.0)
    ap = ActivePowerPI(Kp_p = kp_p, Ki_p = ki_p, ωz = ωz_p, P_ref = p_ref0, ext = Dict{String,Any}("Kω" => kω_droop))
    rq = ReactivePowerPI(Kp_q = kp_q, Ki_q = ki_q, ωf = ωf_q, Q_ref = q_ref0)
    outer = OuterControl(ap, rq)
    inner = CurrentModeControl(kpc = kpc, kic = kic, kffv = kffv)
    dc  = FixedDCSource(voltage = 1200.0)
    pll = KauraPLL(ω_lp = pll_ωlp, kp_pll = pll_kp, ki_pll = pll_ki)
    filt = LCLFilter(lf = lf, rf = rf, cf = cf, lg = lg, rg = rg)
    dyn = DynamicInverter(name = PSY.get_name(g), ω_ref = 1.0, converter = converter,
        outer_control = outer, inner_control = inner, dc_source = dc, freq_estimator = pll, filter = filt)
    add_component_clean!(sys, dyn, g)
    return dyn
end

#Signal extraction
function first_state(results, dev::String, syms::Vector{Symbol})
    for s in syms
        try
            t, x = get_state_series(results, (dev, s))
            return (t, x, s)
        catch
        end
    end
    return nothing
end

normalize_ω(w::AbstractVector) = median(abs.(w)) > 10.0 ? (w ./ Ω_BASE_RAD_S) : w

function dyn_inverter_at_bus(sys::System, busnum::Int)
    target = PSY.get_name(gen_at_bus(sys, busnum))
    for d in get_components(DynamicInverter, sys)
        PSY.get_name(d) == target && return d
    end
    error("No DynamicInverter found at bus $busnum (expected $target)")
end

function get_gfl_pll_speed_pu(results, dev_name::String; kp_pll::Float64 = GFL_PLL_KP, ki_pll::Float64 = GFL_PLL_KI)
    direct = first_state(results, dev_name, Symbol[:ω_pll, :w_pll, :pll_ω, :pll_w, :ω_lp, :w_lp])
    if direct !== nothing
        t, w, src = direct
        return (t, normalize_ω(w), src)
    end
    vq  = first_state(results, dev_name, Symbol[:vq_pll, :v_pll_q, :v_q_pll, :pll_vq])
    eps = first_state(results, dev_name, Symbol[:ε_pll, :eps_pll, :ϵ_pll, :pll_eps])
    if vq !== nothing && eps !== nothing
        t, vqv, _ = vq
        _, epsv, _ = eps
        return (t, 1.0 .+ (kp_pll .* vqv .+ ki_pll .* epsv), :ω_pll_reconstructed)
    end
    return nothing
end

function frequency_dev_hz(results, sys::System, bus::Int)
    name = PSY.get_name(gen_at_bus(sys, bus))
    if bus in SG_BUSES
        try
            t, w = get_state_series(results, (name, :ω))
            return (t, (normalize_ω(w) .- 1.0) .* F_BASE, "Bus $bus (SG)", :ω)
        catch
        end
    end
    if bus in GFM_BUSES
        s = first_state(results, name, Symbol[:ω_oc, :w_oc, :ω, :w])
        s !== nothing && return (s[1], (normalize_ω(s[2]) .- 1.0) .* F_BASE, "Bus $bus (GFM)", s[3])
    end
    if bus in GFL_BUSES
        pll = get_gfl_pll_speed_pu(results, name)
        pll !== nothing && return (pll[1], (pll[2] .- 1.0) .* F_BASE, "Bus $bus (GFL PLL)", pll[3])
    end
    s = first_state(results, name, Symbol[:ω, :w, :ω_oc, :w_oc])
    s !== nothing && return (s[1], (normalize_ω(s[2]) .- 1.0) .* F_BASE, "Bus $bus", s[3])
    return nothing
end

function _unwrap(a::AbstractVector)
    out = copy(a)
    length(out) <= 1 && return out
    for k in 2:length(out)
        d = out[k] - out[k - 1]
        if d > π
            out[k:end] .-= 2π
        elseif d < -π
            out[k:end] .+= 2π
        end
    end
    return out
end

function first_finite(v::AbstractVector)
    for x in v
        isfinite(x) && return x
    end
    return 0.0
end

#Resampling and interpolation
function uniform_grid(Tf::Float64, dt::Float64, events::Vector{Float64})
    tg = collect(0.0:dt:Tf)
    for tev in events
        any(abs.(tg .- tev) .< 1e-12) || push!(tg, tev)
    end
    sort!(tg)
    return tg
end

function interp1(t::AbstractVector, y::AbstractVector, tq::AbstractVector)
    out = fill(NaN, length(tq))
    length(t) == 0 && return out
    ord = sortperm(t); tt = t[ord]; yy = y[ord]
    for (k, x) in pairs(tq)
        (x < tt[1] || x > tt[end]) && continue
        i = searchsortedlast(tt, x)
        if i >= length(tt)
            out[k] = yy[end]
        elseif isapprox(tt[i], x; atol = 1e-12, rtol = 1e-12)
            out[k] = yy[i]
        else
            α = (x - tt[i]) / (tt[i + 1] - tt[i])
            out[k] = (1.0 - α) * yy[i] + α * yy[i + 1]
        end
    end
    return out
end

maybe_resample(t, y, events) = (g = uniform_grid(TF_SIM, EXPORT_DT, events); (g, interp1(t, y, g)))

# Voltage magnitude extraction
_as_float_vector(y) = eltype(y) <: Complex ? Float64.(abs.(y)) : Float64.(y)

function _find_callable(fn_sym::Symbol)
    isdefined(Main, fn_sym) && return getfield(Main, fn_sym)
    isdefined(PSID, fn_sym) && return getfield(PSID, fn_sym)
    return nothing
end

function _maybe_extract_time_value(out)
    if out isa Tuple && length(out) >= 2
        return Float64.(out[1]), _as_float_vector(out[2])
    end
    if out isa DataFrame
        nms = names(out)
        tcols = ["t_s", "time", "t", "Time"]; ycols = ["V_pu", "Vm_pu", "Vmag_pu", "voltage_pu", "magnitude", "Vm"]
        tc = findfirst(c -> c in nms, tcols)
        yc = findfirst(c -> c in nms, ycols)
        if tc !== nothing && yc !== nothing
            return Float64.(out[:, tcols[tc]]), _as_float_vector(out[:, ycols[yc]])
        end
    end
    return nothing
end

function voltage_magnitude_pu(results, sys::System, bus::Int)
    b = bus_component(sys, bus)
    bus_name = PSY.get_name(b)
    fn_candidates = Symbol[:get_voltage_magnitude_series, :get_bus_voltage_magnitude_series,
                           :get_voltage_series, :get_bus_voltage_series]
    arg_candidates = Any[bus, bus_name, b, (bus,), (bus_name,)]
    for fn_sym in fn_candidates
        fn = _find_callable(fn_sym)
        fn === nothing && continue
        for arg in arg_candidates
            try
                tv = _maybe_extract_time_value(fn(results, arg))
                tv !== nothing && return (tv[1], tv[2], @sprintf("Bus %d voltage from %s", bus, string(fn_sym)), fn_sym)
            catch
            end
        end
    end
    dev_name = PSY.get_name(gen_at_bus(sys, bus))
    s = first_state(results, dev_name, Symbol[:V_t, :Vt, :v_t, :vt, :V_mag, :Vmag, :Vm, :V, :v])
    if s !== nothing
        return (s[1], _as_float_vector(s[2]), @sprintf("Bus %d voltage from state %s", bus, string(s[3])), s[3])
    end
    vd = first_state(results, dev_name, Symbol[:vd, :v_d, :Vd, :V_d, :vd_filter, :vr_filter, :vr_cnv])
    vq = first_state(results, dev_name, Symbol[:vq, :v_q, :Vq, :V_q, :vq_filter, :vi_filter, :vi_cnv])
    if vd !== nothing && vq !== nothing
        t_vd, vd_raw, vd_sym = vd
        t_vq, vq_raw, vq_sym = vq
        vq_on_vd = interp1(t_vq, Float64.(vq_raw), t_vd)
        V = sqrt.(Float64.(vd_raw) .^ 2 .+ vq_on_vd .^ 2)
        return (Float64.(t_vd), V, @sprintf("Bus %d voltage from %s/%s", bus, string(vd_sym), string(vq_sym)),
                Symbol("Vmag_from_" * string(vd_sym) * "_" * string(vq_sym)))
    end
    @warn "No bus voltage magnitude series found" bus = bus dev = dev_name
    return nothing
end

# Angle extraction
function angle_symbol_candidates_for_bus(bus::Int)
    if bus in SG_BUSES
        return Symbol[:δ, :delta, :θ, :theta]
    elseif bus in GFM_BUSES
        return Symbol[:θ_oc, :theta_oc, :θ_olc, :θ, :theta]
    elseif bus in GFL_BUSES
        return Symbol[:θ_pll, :theta_pll, :pll_θ, :θ, :theta]
    end
    return Symbol[:δ, :delta, :θ_oc, :θ_pll, :θ, :theta]
end

function raw_device_angle_rad(results, sys::System, bus::Int)
    name = PSY.get_name(gen_at_bus(sys, bus))
    s = first_state(results, name, angle_symbol_candidates_for_bus(bus))
    if s === nothing
        @warn "No angle state found" bus = bus model = model_label(bus) dev = name
        return nothing
    end
    t, ang_raw, sym = s
    return (t, _unwrap(ang_raw), @sprintf("Bus %d raw angle from %s", bus, string(sym)), sym)
end

function rotor_angle_dev_rad(results, sys::System, bus::Int)
    if bus == SLACK_BUS_NUM
        t = collect(0.0:EXPORT_DT:TF_SIM)
        return (t, zeros(length(t)), "Bus $bus fixed reference = 0", :fixed_network_reference_zero)
    end
    raw = raw_device_angle_rad(results, sys, bus)
    raw === nothing && return nothing
    t, θ, _, sym = raw
    return (t, θ .- first_finite(θ), @sprintf("Bus %d angle relative to Bus %d", bus, SLACK_BUS_NUM),
            Symbol("theta$(bus)_minus_bus$(SLACK_BUS_NUM)"))
end

function print_angle_state_diagnostics(results, sys::System)
    println("\n" * "="^100)
    println("ANGLE STATE DIAGNOSTICS (reference = fixed Bus $SLACK_BUS_NUM)")
    println("="^100)
    @printf("%-8s %-12s %-24s %-24s %-18s %-18s\n", "Bus", "Model", "Device", "State", "Initial [rad]", "Final dev [rad]")
    println("-"^100)
    for b in 1:3
        name = PSY.get_name(gen_at_bus(sys, b))
        if b == SLACK_BUS_NUM
            @printf("%-8d %-12s %-24s %-24s %-18.8e %-18.8e\n", b, model_label(b), name, "fixed_ref", 0.0, 0.0)
            continue
        end
        raw = raw_device_angle_rad(results, sys, b)
        if raw === nothing
            @printf("%-8d %-12s %-24s %-24s %-18s %-18s\n", b, model_label(b), name, "NONE", "NA", "NA")
            continue
        end
        t, θ, _, sym = raw
        θ0 = first_finite(θ)
        @printf("%-8d %-12s %-24s %-24s %-18.8e %-18.8e\n", b, model_label(b), name, string(sym), θ0, θ[end] - θ0)
    end
    println("="^100 * "\n")
end

#Power deviations
function mech_power_dev(results, sys::System, bus::Int)
    g = gen_at_bus(sys, bus); name = PSY.get_name(g)
    ratio = safe_gen_base(g) / SYS_BASE_MVA
    for sym in (:x_g1, :x_g2, :x_g3)
        try
            t, x = get_state_series(results, (name, sym))
            Pm = x .* ratio; Pm0 = first_finite(Pm)
            return (t, Pm .- Pm0, Pm0, safe_gen_base(g), sym, @sprintf("Bus %d (SG) ΔPm", bus))
        catch
        end
    end
    return nothing
end

function gfm_power_dev(results, sys::System, bus::Int)
    name = PSY.get_name(gen_at_bus(sys, bus))
    t = nothing; P = nothing; src = :get_activepower_series
    try
        t_s, P_s = get_activepower_series(results, name); t = Float64.(t_s); P = Float64.(P_s)
    catch
    end
    if t === nothing
        for sym in (:Pel, :P_oc, :Pe, :p_el, :active_power)
            try
                t_s, x_s = get_state_series(results, (name, sym)); t = Float64.(t_s); P = Float64.(x_s); src = sym; break
            catch
            end
        end
    end
    t === nothing && (@warn "No GFM active power series" bus = bus; return nothing)
    P0 = first_finite(P)
    return (t, P .- P0, P0, src, @sprintf("Bus %d (GFM) ΔP [P0=%.6f pu]", bus, P0))
end

function gfl_power_loop(results, sys::System, bus::Int)
    dyn = dyn_inverter_at_bus(sys, bus); name = PSY.get_name(dyn)
    pll = get_gfl_pll_speed_pu(results, name)
    pll === nothing && (@warn "No GFL PLL speed series" bus = bus; return nothing)
    t, ω_pll, ωsym = pll
    ap = PSY.get_active_power_control(PSY.get_outer_control(dyn))
    p_ref = try PSY.get_P_ref(ap) catch; 0.0 end
    Kω = try get(PSY.get_ext(ap), "Kω", 0.0) catch; 0.0 end

#GFL reference power back to system base.
    sys_scale = safe_gen_base(gen_at_bus(sys, bus)) / SYS_BASE_MVA
    p_ref_eff = (p_ref .- Kω .* (ω_pll .- 1.0)) .* sys_scale
    p_oc_series = first_state(results, name, Symbol[:p_oc, :P_oc, :poc, :Poc, :active_power_filter])
    if p_oc_series === nothing
        try
            t_p, p_meas = get_activepower_series(results, name)
            p_oc_series = (Float64.(t_p), Float64.(p_meas), :get_activepower_series)
        catch
        end
    end
    p_oc_series === nothing && (@warn "No GFL p_oc series" bus = bus; return nothing)
    t_poc, p_oc_raw, psrc = p_oc_series
    p_oc = interp1(t_poc, p_oc_raw, t)
    p_err = p_ref_eff .- p_oc
    return (t = t, p_oc = p_oc, Δp_oc = p_oc .- first_finite(p_oc),
            p_ref_eff = p_ref_eff, Δp_ref_eff = p_ref_eff .- first_finite(p_ref_eff),
            p_err = p_err, Δp_err = p_err .- first_finite(p_err),
            ω_pll = ω_pll, Kω = Kω, ωsym = ωsym, psrc = psrc)
end

# Terminal P/Q
function terminal_power_series(results, sys::System, bus::Int; channel::Symbol)
    name = PSY.get_name(gen_at_bus(sys, bus))
    if channel == :P
        try
            t, p = get_activepower_series(results, name); return (Float64.(t), Float64.(p), :get_activepower_series)
        catch
        end
        for sym in (:Pel, :P_el, :Pe, :p_el, :active_power, :P, :p)
            try
                t, x = get_state_series(results, (name, sym)); return (Float64.(t), Float64.(x), sym)
            catch
            end
        end
        return nothing
    else
        try
            t, q = get_reactivepower_series(results, name); return (Float64.(t), Float64.(q), :get_reactivepower_series)
        catch
        end
        for sym in (:Qel, :Q_el, :Qe, :q_el, :reactive_power, :Q, :q)
            try
                t, x = get_state_series(results, (name, sym)); return (Float64.(t), Float64.(x), sym)
            catch
            end
        end
        return nothing
    end
end

function terminal_pq_dev(results, sys::System, bus::Int)
    ps = terminal_power_series(results, sys, bus; channel = :P)
    qs = terminal_power_series(results, sys, bus; channel = :Q)
    (ps === nothing && qs === nothing) && return nothing
    if ps !== nothing
        tP, P, psrc = ps; t_ref = tP
    else
        t_ref = qs[1]; P = fill(NaN, length(t_ref)); psrc = :missing
    end
    if qs !== nothing
        tQ, Q, qsrc = qs
    else
        Q = fill(NaN, length(t_ref)); qsrc = :missing
    end
    (ps !== nothing && qs !== nothing) && (Q = interp1(tQ, Q, t_ref))
    P0 = first_finite(P); Q0 = first_finite(Q)
    return (t = t_ref, P = P, Q = Q, ΔP = P .- P0, ΔQ = Q .- Q0, P0 = P0, Q0 = Q0,
            psrc = psrc, qsrc = qsrc, dev = PSY.get_name(gen_at_bus(sys, bus)))
end

function save_terminal_pq_csv(results, sys, bus, outdir, events)
    mkpath(outdir)
    out = terminal_pq_dev(results, sys, bus)
    out === nothing && (@warn "No terminal P/Q for CSV" bus = bus; return nothing)
    t, ΔP   = maybe_resample(out.t, out.ΔP, events)
    _, ΔQ   = maybe_resample(out.t, out.ΔQ, events)
    _, Pabs = maybe_resample(out.t, out.P, events)
    _, Qabs = maybe_resample(out.t, out.Q, events)
    path = joinpath(outdir, @sprintf("terminal_pq_bus%02d.csv", bus))
    CSV.write(path, DataFrame(t_s = t, delta_p_pu = ΔP, p_abs_pu = Pabs, p_initial_pu = fill(out.P0, length(t)),
        delta_q_pu = ΔQ, q_abs_pu = Qabs, q_initial_pu = fill(out.Q0, length(t)),
        asset_type = fill(model_label(bus), length(t)), device_name = fill(out.dev, length(t))))
    return path
end

function print_terminal_pq_report(results, sys; tail_s = 2.0)
    println("\n" * "="^112)
    println("TERMINAL P/Q INITIALIZATION AND FINAL-TAIL CHECK")
    println("="^112)
    @printf("%-6s %-10s %-14s %-14s %-14s %-14s %-14s %-14s\n",
            "Bus", "Model", "P0 [pu]", "Q0 [pu]", "ΔP_last", "ΔQ_last", "ΔP_tail", "ΔQ_tail")
    println("-"^112)
    for bus in sort(unique([SG_BUSES; GFM_BUSES; GFL_BUSES]))
        out = terminal_pq_dev(results, sys, bus)
        if out === nothing
            @printf("%-6d %-10s %s\n", bus, model_label(bus), "(no terminal P/Q)"); continue
        end
        i0 = findfirst(>=(out.t[end] - tail_s), out.t); i0 = i0 === nothing ? 1 : i0
        @printf("%-6d %-10s %-14.6f %-14.6f %-14.6f %-14.6f %-14.6f %-14.6f\n",
                bus, model_label(bus), out.P0, out.Q0, out.ΔP[end], out.ΔQ[end],
                mean(skipmissing(out.ΔP[i0:end])), mean(skipmissing(out.ΔQ[i0:end])))
    end
    println("="^112 * "\n")
end

#Network P/Q calculation
function stamp_branch_ybus!(Y::AbstractMatrix{ComplexF64}, fbus::Int, tbus::Int,
                            r::Float64, x::Float64, b::Float64, tap_raw::Float64;
                            buses::Vector{Int})
    bus_to_idx = Dict{Int,Int}(bus => i for (i, bus) in enumerate(buses))
    (haskey(bus_to_idx, fbus) && haskey(bus_to_idx, tbus)) || return Y

    z = complex(r, x)
    abs(z) < 1e-12 && return Y

    tap = abs(tap_raw) < 1e-12 ? 1.0 : tap_raw
    y = inv(z)
    ysh = 1im * b / 2.0
    i = bus_to_idx[fbus]
    j = bus_to_idx[tbus]

    Y[i, i] += (y + ysh) / (tap * tap)
    Y[i, j] -= y / tap
    Y[j, i] -= y / tap
    Y[j, j] += y + ysh
    return Y
end

function build_ybus_from_matpower(; buses::Vector{Int} = collect(1:3))
    Y = zeros(ComplexF64, length(buses), length(buses))
    for (fbus, tbus, r, x, b, tap) in bus3_BRANCH_DATA
        stamp_branch_ybus!(Y, fbus, tbus, r, x, b, tap; buses = buses)
    end
    return Y
end

function build_ybus_dense_from_system(sys::System; buses::Vector{Int} = collect(1:3))
    bus_to_idx = Dict{Int,Int}(bus => i for (i, bus) in enumerate(buses))
    Y = zeros(ComplexF64, length(buses), length(buses))
    nstamped = 0

    for line in get_components(Line, sys)
        fbus, tbus = line_pair(line)
        (haskey(bus_to_idx, fbus) && haskey(bus_to_idx, tbus)) || continue

        r = try Float64(PSY.get_r(line)) catch; NaN end
        x = try Float64(PSY.get_x(line)) catch; NaN end
        b = try Float64(PSY.get_b(line)) catch; 0.0 end
        (!isfinite(r) || !isfinite(x)) && continue

        stamp_branch_ybus!(Y, fbus, tbus, r, x, b, 0.0; buses = buses)
        nstamped += 1
    end

    if nstamped == 0
        @info "Building Ybus from MATPOWER branch table"
        return build_ybus_from_matpower(; buses = buses)
    end

    return Y
end

function pf_reference_vectors(pf::PowerFlowData; buses::Vector{Int} = collect(1:3))
    V0 = zeros(Float64, length(buses)); θ0 = zeros(Float64, length(buses))
    for (i, b) in enumerate(buses)
        idx = pf.bus_lookup[b]; V0[i] = Float64(pf.bus_magnitude[idx]); θ0[i] = Float64(pf.bus_angles[idx])
    end
    θ0 .-= θ0[end]
    return V0, θ0
end

function local_load_step_series(t::AbstractVector, steps::Vector{RealizedLoadStep}, bus::Int)
    ΔP = zeros(Float64, length(t)); ΔQ = zeros(Float64, length(t))
    for s in steps
        s.bus == bus || continue
        for k in eachindex(t)
            if t[k] >= s.t_step - 1e-12
                ΔP[k] += s.deltaP_pu; ΔQ[k] += s.deltaQ_pu
            end
        end
    end
    return ΔP, ΔQ
end

function network_bus_pq_dev(results, sys::System, pf::PowerFlowData, steps::Vector{RealizedLoadStep})
    buses = collect(1:3)
    t_grid = uniform_grid(TF_SIM, EXPORT_DT, unique(sort([s.t_step for s in steps])))
    Vmag = zeros(Float64, length(buses), length(t_grid))
    θdev = zeros(Float64, length(buses), length(t_grid))
    angle_source = String[]
    for (i, b) in enumerate(buses)
        vout = voltage_magnitude_pu(results, sys, b)
        vout === nothing && (@warn "Missing voltage for network P/Q" bus = b; return nothing)
        tv, V, _, _ = vout
        Vmag[i, :] .= interp1(tv, V, t_grid)
        if b == SLACK_BUS_NUM
            push!(angle_source, "fixed_reference_zero"); continue
        end
        aout = raw_device_angle_rad(results, sys, b)
        aout === nothing && (@warn "Missing angle for network P/Q" bus = b; return nothing)
        ta, θ, _, θsym = aout
        θg = interp1(ta, θ, t_grid)
        θdev[i, :] .= θg .- first_finite(θg)
        push!(angle_source, "device_angle_" * string(θsym))
    end
    Y = build_ybus_dense_from_system(sys; buses = buses)
    V0, θ0 = pf_reference_vectors(pf; buses = buses)
    Sref = (V0 .* exp.(1im .* θ0)) .* conj.(Y * (V0 .* exp.(1im .* θ0)))
    P = zeros(Float64, length(buses), length(t_grid)); Q = zeros(Float64, length(buses), length(t_grid))
    for k in eachindex(t_grid)
        Vc = Vmag[:, k] .* exp.(1im .* (θ0 .+ θdev[:, k]))
        S = Vc .* conj.(Y * Vc)
        P[:, k] .= real.(S); Q[:, k] .= imag.(S)
    end
    @warn "Network P/Q is diagnostic (angle reconstruction; excludes shunt loads)." angle_source = angle_source
    return (t = t_grid, buses = buses, P = P, Q = Q, Pref = real.(Sref), Qref = imag.(Sref))
end

function save_network_pq_csv(results, sys::System, pf::PowerFlowData, bus::Int, outdir::String, steps::Vector{RealizedLoadStep})
    mkpath(outdir)
    net = network_bus_pq_dev(results, sys, pf, steps)
    net === nothing && (@warn "No network P/Q for CSV" bus = bus; return nothing)
    idx = findfirst(==(bus), net.buses); idx === nothing && return nothing
    Pabs = vec(net.P[idx, :]); Qabs = vec(net.Q[idx, :])
    ΔP_net = Pabs .- net.Pref[idx]; ΔQ_net = Qabs .- net.Qref[idx]
    ΔP_load, ΔQ_load = local_load_step_series(net.t, steps, bus)
    path = joinpath(outdir, @sprintf("network_pq_bus%02d.csv", bus))
    CSV.write(path, DataFrame(t_s = net.t, delta_p_pu = ΔP_net .+ ΔP_load, delta_q_pu = ΔQ_net .+ ΔQ_load,
        delta_p_network_pu = ΔP_net, delta_q_network_pu = ΔQ_net,
        p_network_abs_pu = Pabs, q_network_abs_pu = Qabs,
        delta_p_load_pu = ΔP_load, delta_q_load_pu = ΔQ_load,
        asset_type = fill(model_label(bus), length(net.t))))
    return path
end

# CSV exports
function save_freq_csv(results, sys, bus, outdir, events)
    mkpath(outdir)
    s = frequency_dev_hz(results, sys, bus)
    s === nothing && (@warn "No frequency series" bus = bus; return nothing)
    t, df = maybe_resample(s[1], s[2], events)
    path = joinpath(outdir, @sprintf("freq_bus%02d.csv", bus))
    CSV.write(path, DataFrame(t_s = t, df_hz = df)); return path
end

function save_angle_csv(results, sys, bus, outdir, events)
    mkpath(outdir)
    s = rotor_angle_dev_rad(results, sys, bus)
    s === nothing && (@warn "No angle series" bus = bus; return nothing)
    t, dδ = maybe_resample(s[1], s[2], events)
    path = joinpath(outdir, @sprintf("rotor_angle_bus%02d.csv", bus))
    CSV.write(path, DataFrame(t_s = t, delta_rad = dδ)); return path
end

function save_voltage_mag_csv(results, sys, bus, outdir, events)
    mkpath(outdir)
    out = voltage_magnitude_pu(results, sys, bus)
    out === nothing && (@warn "No voltage series" bus = bus; return nothing)
    t, V = maybe_resample(out[1], out[2], events); V0 = first_finite(V)
    path = joinpath(outdir, @sprintf("voltage_bus%02d.csv", bus))
    CSV.write(path, DataFrame(t_s = t, V_pu = V, delta_V_pu = V .- V0, V_initial_pu = fill(V0, length(t)))); return path
end

function save_pg_csv(results, sys, bus, outdir, events)
    mkpath(outdir)
    out = mech_power_dev(results, sys, bus)
    out === nothing && (@warn "No TGTypeI mechanical power" bus = bus; return nothing)
    t, ΔPm, Pm0, _, _, _ = out
    t, ΔPm = maybe_resample(t, ΔPm, events)
    path = joinpath(outdir, @sprintf("pg_bus%02d.csv", bus))
    CSV.write(path, DataFrame(t_s = t, delta_pm_pu = ΔPm, Pm0_pu_sysbase = fill(Pm0, length(t)))); return path
end

function save_gfm_power_csv(results, sys, bus, outdir, events)
    mkpath(outdir)
    out = gfm_power_dev(results, sys, bus)
    out === nothing && (@warn "No GFM power series" bus = bus; return nothing)
    t, ΔP, P0, _, _ = out
    t, ΔP = maybe_resample(t, ΔP, events)
    path = joinpath(outdir, @sprintf("gfm_power_dev_bus%02d.csv", bus))
    CSV.write(path, DataFrame(t_s = t, delta_p_pu = ΔP, p_abs_pu = ΔP .+ P0, p_initial_pu = fill(P0, length(t)))); return path
end

function save_gfl_power_loop_csv(results, sys, bus, outdir, events)
    mkpath(outdir)
    out = gfl_power_loop(results, sys, bus)
    out === nothing && (@warn "No GFL power loop series" bus = bus; return nothing)
    t, p_oc       = maybe_resample(out.t, out.p_oc, events)
    _, Δp_oc      = maybe_resample(out.t, out.Δp_oc, events)
    _, p_ref_eff  = maybe_resample(out.t, out.p_ref_eff, events)
    _, p_err      = maybe_resample(out.t, out.p_err, events)
    _, ω_pll      = maybe_resample(out.t, out.ω_pll, events)
    path = joinpath(outdir, @sprintf("gfl_power_loop_bus%02d.csv", bus))
    CSV.write(path, DataFrame(t_s = t, p_oc_pu = p_oc, delta_p_oc_pu = Δp_oc, p_ref_eff_pu = p_ref_eff,
        p_error_pu = p_err, omega_pll_pu = ω_pll, Kw = fill(out.Kω, length(t)))); return path
end

# Power-flow export
function export_pf_operating_point(pf::PowerFlowData, outdir::String)
    mkpath(outdir)
    println("\n--- Power-flow operating point ---")
    @printf("%-6s %-10s %-10s %-12s %-12s\n", "Bus", "Vm [pu]", "Va [deg]", "P_net [MW]", "Q_net [MVAr]")
    println("-"^62)
    rows = DataFrame(bus = Int[], Vm_pu = Float64[], Va_deg = Float64[], P_net_MW = Float64[], Q_net_MVAr = Float64[])
    for b in 1:3
        haskey(pf.bus_lookup, b) || continue
        idx = pf.bus_lookup[b]
        Vm = pf.bus_magnitude[idx]; Va = pf.bus_angles[idx] * 180.0 / π
        Pnet = (safe_pinj(pf, idx) - safe_pwd(pf, idx)) * SYS_BASE_MVA
        Qnet = (safe_qinj(pf, idx) - safe_qwd(pf, idx)) * SYS_BASE_MVA
        @printf("%-6d %-10.4f %-10.4f %-12.4f %-12.4f\n", b, Vm, Va, Pnet, Qnet)
        push!(rows, (b, Vm, Va, Pnet, Qnet))
    end
    CSV.write(joinpath(outdir, "powerflow_operating_point.csv"), rows)
    return rows
end

# Steady-state frequency report
function print_ss_freq(results, sys; tail_s = 2.0)
    println("\n" * "="^80)
    println("STEADY-STATE FREQUENCY DEVIATION (last $(tail_s) s)")
    println("="^80)
    @printf("%-6s %-12s %-14s %-14s %-14s %-8s\n", "Bus", "Model", "Δf_last[Hz]", "Δf_mean[Hz]", "Δf_std[Hz]", "Ntail")
    println("-"^80)
    for bus in 1:3
        s = frequency_dev_hz(results, sys, bus)
        if s === nothing
            @printf("%-6d %-12s %s\n", bus, model_label(bus), "(no series)"); continue
        end
        t, df, lab, _ = s
        i0 = findfirst(>=(t[end] - tail_s), t); i0 = i0 === nothing ? 1 : i0
        tail = df[i0:end]
        @printf("%-6d %-12s %-14.6f %-14.6f %-14.6f %-8d\n", bus, model_label(bus), df[end], mean(tail), std(tail), length(tail))
    end
    println("="^80 * "\n")
end

# Small-signal analysis
function extract_eigenvalues(small_sig)
    for fld in (:eigenvalues, :eigvals, :λ, :lambda)
        try return ComplexF64.(collect(getproperty(small_sig, fld))) catch end
    end
    error("Could not extract eigenvalues from small_signal_analysis result.")
end

function damping_pct(λ::ComplexF64)
    mag = abs(λ)
    mag <= eps(Float64) && return 100.0
    return max(0.0, -real(λ) / mag * 100.0)
end

classify_mode(λ::ComplexF64) = real(λ) > SSA_STABILITY_TOL ? "Unstable" :
                               abs(real(λ)) <= SSA_REAL_NEAR_ZERO ? "Near-axis" : "Stable"
mode_type(λ::ComplexF64) = abs(imag(λ)) > 1e-4 ? "Osc" : "NonOsc"

function make_ssa_mode_dataframe(eigs::Vector{ComplexF64})
    rows = NamedTuple[]
    for (i, λ) in enumerate(eigs)
        ζ = damping_pct(λ); status = classify_mode(λ); typ = mode_type(λ)
        push!(rows, (mode = i, real_part = real(λ), imag_part = imag(λ), abs_value = abs(λ),
                     freq_hz = abs(imag(λ)) / (2π), damping_pct = ζ, status = status, type = typ,
                     unstable = status == "Unstable", near_axis = status == "Near-axis",
                     lightly_damped_osc = typ == "Osc" && ζ < SSA_DAMPING_WARN_PCT))
    end
    return DataFrame(rows)
end

function plot_ssa_eigs(mode_df::DataFrame, fig_path::String)
    p = plot(xlabel = "Real(λ)", ylabel = "Imag(λ)", title = "SSA Eigenvalues",
             legend = :best, grid = false, framestyle = :box, dpi = 150, size = (950, 650),
             titlefont = font(12, :bold), guidefont = font(13, :bold), tickfont = font(11, :bold),
             xlims = (-4000, 4000), ylims = (-4000, 4000))
    sdf = mode_df[mode_df.unstable .== false, :]
    udf = mode_df[mode_df.unstable .== true, :]
    nrow(sdf) > 0 && scatter!(p, sdf.real_part, sdf.imag_part; ms = 4, label = "Stable / near-axis")
    nrow(udf) > 0 && scatter!(p, udf.real_part, udf.imag_part; ms = 6, marker = :xcross, label = "Unstable")
    vline!(p, [0.0]; lw = 2, ls = :dash, c = :black, alpha = 0.75, label = "Re(λ)=0")
    savefig(p, fig_path)
    return fig_path
end

function write_ssa_txt_report(mode_df::DataFrame, path::String; label::String = SSA_REPORT_LABEL)
    open(path, "w") do io
        println(io, "="^110)
        println(io, "SSA STABILITY REPORT  |  $(Dates.now())  |  $label")
        println(io, "Total modes: $(nrow(mode_df))")
        println(io, "="^110, "\n")
        rightmost = mode_df[argmax(mode_df.real_part), :]
        @printf(io, "Rightmost eigenvalue : %.8e %+.8ej  (mode %d)\n", rightmost.real_part, rightmost.imag_part, rightmost.mode)
        @printf(io, "Unstable modes       : %d\n", sum(mode_df.unstable))
        @printf(io, "Near-axis modes      : %d\n", sum(mode_df.near_axis))
        @printf(io, "Lightly damped osc.  : %d\n\n", sum(mode_df.lightly_damped_osc))
        @printf(io, "%-8s %-26s %-12s %-12s %-12s %-8s\n", "Mode", "Eigenvalue", "Freq[Hz]", "Damp[%]", "Status", "Type")
        println(io, "-"^110)
        for r in eachrow(sort(mode_df, [:real_part], rev = true))
            @printf(io, "%-8d %-26s %-12.5f %-12.5f %-12s %-8s\n",
                    r.mode, @sprintf("%.6e %+.6ej", r.real_part, r.imag_part), r.freq_hz, r.damping_pct, r.status, r.type)
        end
        println(io, "="^110)
    end
    return path
end

function run_ssa_reporting(sim, outdir::String; label::String = SSA_REPORT_LABEL)
    mkpath(outdir)
    small_sig = with_logger(NullLogger()) do
        small_signal_analysis(sim)
    end
    eigs = extract_eigenvalues(small_sig)
    mode_df = make_ssa_mode_dataframe(eigs)
    CSV.write(joinpath(outdir, SSA_EIG_CSV_FILENAME), mode_df)
    CSV.write(joinpath(outdir, SSA_SUMMARY_CSV_FILENAME), sort(mode_df, [:real_part], rev = true))
    write_ssa_txt_report(mode_df, joinpath(outdir, SSA_TXT_FILENAME); label = label)
    write_ssa_txt_report(mode_df, joinpath(outdir, SSA_FULL_MODAL_REPORT_FILENAME); label = label)
    plot_ssa_eigs(mode_df, joinpath(outdir, SSA_PLOT_FILENAME))
    rightmost = mode_df[argmax(mode_df.real_part), :]
    println("\n" * "="^100)
    println("SSA SUMMARY")
    println("="^100)
    @printf("Rightmost mode: %d | λ = %.8e %+.8ej | f = %.5f Hz | ζ = %.5f %% | %s\n",
            rightmost.mode, rightmost.real_part, rightmost.imag_part, rightmost.freq_hz, rightmost.damping_pct, rightmost.status)
    println("Unstable modes: $(sum(mode_df.unstable))")
    println("SSA outputs saved to: $outdir")
    println("="^100 * "\n")
    return (small_sig = small_sig, eigs = eigs, mode_df = mode_df)
end

# Solver

_status_text(sim, ret) =
    uppercase(string(try string(ret) catch; "" end, " ",
                     try string(getproperty(sim, :status)) catch; "" end))

function solution_final_time(results, sys::System)
    for bus in 1:3
        s = try frequency_dev_hz(results, sys, bus) catch; nothing end
        (s !== nothing && !isempty(s[1])) && return Float64(s[1][end])
    end
    return -Inf
end

function execute_with_fallback!(build_sim::Function, sys::System; tfinal::Float64 = TF_SIM)
    attempts = [
        (name = "IDA(ResidualModel)", mtype = ResidualModel,
         solver = IDA(linear_solver = :Dense, max_order = SOLVER_MAX_ORDER),
         kwargs = (abstol = SOLVER_ABSTOL, reltol = SOLVER_RELTOL, dtmax = SOLVER_DTMAX)),
        (name = "Rodas5P(MassMatrix)", mtype = MassMatrixModel,
         solver = Rodas5P(autodiff = false),
         kwargs = (abstol = SOLVER_ABSTOL, reltol = SOLVER_RELTOL, dt = SOLVER_INIT_DT, dtmax = SOLVER_DTMAX)),
        (name = "FBDF(MassMatrix)", mtype = MassMatrixModel,
         solver = FBDF(autodiff = false),
         kwargs = (abstol = SOLVER_ABSTOL, reltol = SOLVER_RELTOL, dt = SOLVER_INIT_DT, dtmax = SOLVER_DTMAX)),
    ]
    last_err = nothing
    for a in attempts
        sim = build_sim(a.mtype)
        try
            ret = execute!(sim, a.solver; a.kwargs...)
            st = _status_text(sim, ret)
            if occursin("FAIL", st) || occursin("UNSTABLE", st) || occursin("INCOMPLETE", st)
                @warn "$(a.name) reported a non-finalized status; trying next solver." status = st
                continue
            end
            res = read_results(sim)
            tend = solution_final_time(res, sys)
            if res !== nothing && tend >= 0.95 * tfinal
                println("Integrated with $(a.name) (reached t = $(round(tend, digits = 3)) s).")
                return sim, res
            end
            @warn "$(a.name) aborted before the final time; trying next solver." reached_s = tend target_s = tfinal
        catch err
            last_err = err
            @warn "$(a.name) failed; trying next solver." err = string(err)
        end
    end
    error("All solvers failed to reach t = $(tfinal) s. The model initializes and is " *
          "small-signal stable, so suspect the load-step event handling — try a smaller " *
          "step or a constant-impedance perturbation. Last error: $(last_err)")
end

# Plotting
function _new_panel(title, ylabel; xlabel = "")
    return plot(title = title, ylabel = ylabel, xlabel = xlabel, legend = :best,
                grid = false, framestyle = :box, xlims = (0.0, TF_SIM), widen = false,
                dpi = 140, size = (1400, 500), titlefont = font(12, :bold),
                guidefont = font(14, :bold), tickfont = font(12, :bold), legendfont = font(11))
end

function _mark_events!(p, ev_times)
    for tev in ev_times
        vline!(p, [tev]; lw = 1, ls = :dash, c = :black, alpha = 0.35, label = "")
    end
end

function plot_voltage_magnitude_panel(results, sys, steps, plotdir)
    mkpath(plotdir)
    ev = unique(sort([s.t_step for s in steps]))
    p = _new_panel("Bus voltage magnitudes\n" * event_summary(steps), "V [pu]"; xlabel = "Time [s]")
    plot!(p; ylims = ABS_VOLTAGE_YLIM, yticks = ABS_VOLTAGE_YTICKS)
    for bus in sort(unique([SG_BUSES; GFM_BUSES; GFL_BUSES]))
        out = voltage_magnitude_pu(results, sys, bus); out === nothing && continue
        t, V = maybe_resample(out[1], out[2], ev)
        plot!(p, t, V;  lw = 3.0, color = bus_plot_color(bus), label = bus_label(bus))
    end
    _mark_events!(p, ev)
    savefig(p, joinpath(plotdir, "3bus_voltage_magnitude_pu.png"))
    return p
end

function plot_voltage_deviation_panel(results, sys, steps, plotdir)
    mkpath(plotdir)
    ev = unique(sort([s.t_step for s in steps]))
    p = _new_panel("Bus voltage-magnitude deviations\n" * event_summary(steps), "ΔV [pu]"; xlabel = "Time [s]")
    for bus in sort(unique([SG_BUSES; GFM_BUSES; GFL_BUSES]))
        out = voltage_magnitude_pu(results, sys, bus); out === nothing && continue
        t, V = maybe_resample(out[1], out[2], ev); V0 = first_finite(V)
        plot!(p, t, V .- V0; ylims = (-0.06, 0.00), lw = 3.0, color = bus_plot_color(bus), label = @sprintf("%s (V₀=%.4f)", bus_label(bus), V0))
    end
    hline!(p, [0.0]; lw = 0.8, ls = :dot, c = :gray, label = "")
    _mark_events!(p, ev)
    savefig(p, joinpath(plotdir, "3bus_voltage_deviation_pu.png"))
    return p
end

function plot_summary(results, sys, steps, plotdir)
    mkpath(plotdir)
    ev = unique(sort([s.t_step for s in steps]))
    ev_str = event_summary(steps)
    all_buses = sort(unique([SG_BUSES; GFM_BUSES; GFL_BUSES]))

    p1 = _new_panel("Frequency state\n" * ev_str, "Δf [Hz]")
    for bus in all_buses
        s = frequency_dev_hz(results, sys, bus); s === nothing && continue
        t, df = maybe_resample(s[1], s[2], ev)
        plot!(p1, t, df; ylims = (-0.06, 0.00), lw = 3.0, color = bus_plot_color(bus), label = bus_label(bus))
    end
    _mark_events!(p1, ev)

    p2 = _new_panel("Angle states (Δθ from initial)", "Δθ [rad]")
    for bus in all_buses
        s = rotor_angle_dev_rad(results, sys, bus); s === nothing && continue
        t, dδ = maybe_resample(s[1], s[2], ev)
        plot!(p2, t, dδ; ylims = (-0.07, 0.00), lw = 3.0, color = bus_plot_color(bus), label = bus_label(bus))
    end
    _mark_events!(p2, ev)

    p3 = _new_panel("Terminal active power", "ΔP [pu, system base]")
    for bus in all_buses
        out = terminal_pq_dev(results, sys, bus); out === nothing && continue
        t, ΔP = maybe_resample(out.t, out.ΔP, ev)
        plot!(p3, t, ΔP; ylims = (-0.2, 0.20), lw = 3.0, color = bus_plot_color(bus), label = @sprintf("%s (P₀=%.4f)", bus_label(bus), out.P0))
    end
    hline!(p3, [0.0]; lw = 0.8, ls = :dot, c = :gray, label = "")
    _mark_events!(p3, ev)

    p4 = _new_panel("Terminal reactive power", "ΔQ [pu, system base]")
    for bus in all_buses
        out = terminal_pq_dev(results, sys, bus); out === nothing && continue
        t, ΔQ = maybe_resample(out.t, out.ΔQ, ev)
        plot!(p4, t, ΔQ; ylims = (-0.06, 0.08), lw = 3.0, color = bus_plot_color(bus), label = @sprintf("%s (Q₀=%.4f)", bus_label(bus), out.Q0))
    end
    hline!(p4, [0.0]; lw = 0.8, ls = :dot, c = :gray, label = "")
    _mark_events!(p4, ev)

    p5 = _new_panel("SG mechanical power state", "ΔPₘ [pu, system base]")
    for bus in SG_BUSES
        out = mech_power_dev(results, sys, bus); out === nothing && continue
        t, ΔPm, Pm0, mbase, _, _ = out
        t, ΔPm = maybe_resample(t, ΔPm, ev)
        plot!(p5, t, ΔPm; ylims = (0.00, 0.06), lw = 3.0, color = bus_plot_color(bus),
              label = @sprintf("%s (mBase=%.0f, Pm₀=%.3f)", bus_label(bus), mbase, Pm0))
    end
    hline!(p5, [0.0]; lw = 0.8, ls = :dot, c = :gray, label = "")
    _mark_events!(p5, ev)

    p6 = _new_panel("GFL active-power loop: measured p_oc", "Δp_oc [pu]"; xlabel = "Time [s]")
    for bus in GFL_BUSES
        out = gfl_power_loop(results, sys, bus); out === nothing && continue
        t, Δp_oc = maybe_resample(out.t, out.Δp_oc, ev)
        plot!(p6, t, Δp_oc; ylims = (-0.012, 0.015), lw = 3.0, color = bus_plot_color(bus), label = @sprintf("%s (Kω=%.2f)", bus_label(bus), out.Kω))
    end
    hline!(p6, [0.0]; lw = 0.8, ls = :dot, c = :gray, label = "")
    _mark_events!(p6, ev)

    fig = plot(p1, p2, p3, p4, p5, p6; layout = (6, 1), size = (1500, 2300), dpi = 140,
               left_margin = 10Plots.mm, right_margin = 4Plots.mm, bottom_margin = 6Plots.mm)
    savefig(fig, joinpath(plotdir, "3bus_summary.png"))
    savefig(p1, joinpath(plotdir, "3bus_frequency.png"))
    savefig(p2, joinpath(plotdir, "3bus_relative_angle.png"))
    savefig(p3, joinpath(plotdir, "3bus_terminal_active_power.png"))
    savefig(p4, joinpath(plotdir, "3bus_terminal_reactive_power.png"))
    savefig(p5, joinpath(plotdir, "3bus_sg_mechanical_power.png"))
    savefig(p6, joinpath(plotdir, "3bus_gfl_active_power_loop.png"))
    return fig
end

#Main Run
function main(; show_system::Bool = true)
    validate_partitions!(SG_BUSES, GFM_BUSES, GFL_BUSES)

    println("="^80)
    println("3-bus dq0 EMT simulation with SG / GFM / GFL resources")
    println("SG = $SG_BUSES   GFM = $GFM_BUSES   GFL = $GFL_BUSES   slack = $SLACK_BUS_NUM")
    println("="^80)

    steps = LOAD_CHANGE_ENABLE ? build_realized_load_steps() : RealizedLoadStep[]
    isempty(steps) || print_load_schedule(steps)
    ev_times = unique(sort([s.t_step for s in steps]))

    plotdir = joinpath(pwd(), "3bus_plots")
    csvdir  = joinpath(plotdir, "csv")
    pfdir   = joinpath(plotdir, "pf")
    ssadir  = joinpath(plotdir, "ssa")
    foreach(mkpath, (plotdir, csvdir, pfdir, ssadir))

    case_path = joinpath(pwd(), "modified_3bus_configurable.m")
    write_modified_3bus_case(case_path)

    sys = System(case_path; runchecks = false)
    PSY.set_units_base_system!(sys, "SYSTEM_BASE")
    enforce_load_bases!(sys)

    selected_pairs = selected_dynamic_line_pair_set(sys)
    print_dynamic_line_selection(sys, selected_pairs)
    print_dynamic_line_damping_edits(apply_dynamic_line_damping_edits!(sys, selected_pairs))
    add_dynamic_lines!(sys; include_pairs = selected_pairs)

    if show_system
        println("\n--- System summary before dynamic attachment ---")
        display(sys); println()
    end

    pf_data = PowerFlowData(ACPowerFlow(), sys)
    solve_powerflow!(pf_data)
    update_system_voltages!(sys, pf_data)

    loads_to_constant_impedance!(sys)

    sync_static_gens_to_pf!(sys, pf_data; buses = collect(1:3))
    set_generator_bases!(sys)

    export_pf_operating_point(pf_data, pfdir)
    print_reference_summary(sys)

    mach_p = sg_machine_params()
    avr_p  = sg_avr_params()
    for b in SG_BUSES
        attach_sg!(sys, gen_at_bus(sys, b); H = getH(b), D = getDsg(b), droop = R_SG_PU, mach_p = mach_p, avr_p = avr_p)
    end
    for b in GFM_BUSES
        attach_gfm!(sys, gen_at_bus(sys, b); droop = R_GFM_PU, Ta = getGfmTa(b), kd = getGfmKd(b), kω = getGfmKw(b))
    end
    for b in GFL_BUSES
        attach_gfl!(sys, gen_at_bus(sys, b); kω_droop = getGflKw(b))
    end

    for b in 1:3
        ensure_perturb_load!(sys, b)
    end
    perts = LOAD_CHANGE_ENABLE ? make_load_changes_from_event_schedule(sys, steps) : PSID.Perturbation[]

    build_sim(mtype) = Simulation!(mtype, sys, pwd(), (0.0, TF_SIM), perts)

    if RUN_SSA_REPORT
        try
            run_ssa_reporting(build_sim(ResidualModel), ssadir; label = SSA_REPORT_LABEL)
        catch err
            @warn "SSA reporting failed; time-domain simulation continues" err = string(err)
        end
    end

    sim, results = execute_with_fallback!(build_sim, sys)
    results === nothing && (@error "No results returned."; return nothing)

    print_angle_state_diagnostics(results, sys)
    print_ss_freq(results, sys)
    print_terminal_pq_report(results, sys)

    for b in 1:3
        save_freq_csv(results, sys, b, csvdir, ev_times)
        save_angle_csv(results, sys, b, csvdir, ev_times)
        save_voltage_mag_csv(results, sys, b, csvdir, ev_times)
        save_terminal_pq_csv(results, sys, b, csvdir, ev_times)
        SAVE_DIAGNOSTIC_NETWORK_PQ && save_network_pq_csv(results, sys, pf_data, b, csvdir, steps)
    end
    for b in SG_BUSES;  save_pg_csv(results, sys, b, csvdir, ev_times) end
    for b in GFM_BUSES; save_gfm_power_csv(results, sys, b, csvdir, ev_times) end
    for b in GFL_BUSES; save_gfl_power_loop_csv(results, sys, b, csvdir, ev_times) end

    plot_summary(results, sys, steps, plotdir)
    plot_voltage_magnitude_panel(results, sys, steps, plotdir)
    plot_voltage_deviation_panel(results, sys, steps, plotdir)

    println("Done. Plots in $plotdir, CSVs in $csvdir")
    return (sys = sys, results = results, steps = steps, plotdir = plotdir, csvdir = csvdir, ssadir = ssadir)
end

main(show_system = true)
