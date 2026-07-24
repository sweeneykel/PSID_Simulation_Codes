#3bus_parameters.jl

#settings for the 3-bus SG/GFM/GFL study.

const SYS_BASE_MVA = 100.0   #System base[MVA]
const F_BASE       = 60.0     #Nominal frequency[Hz]
const Ω_BASE_RAD_S = 2π * F_BASE     #Nominal speed[rad/s]

const TF_SIM    = 10.0        #Simulation time[s]
const EXPORT_DT = 0.01   #data Export time step[s]

# Solver settings
const SOLVER_ABSTOL    = 1e-6     #IDA absolute tolerance
const SOLVER_RELTOL    = 1e-6         #IDA relative tolerance
const SOLVER_MAX_ORDER = 5      #IDA maximum order
const SOLVER_INIT_DT   = 1e-4     #IDA initial step[s]
const SOLVER_DTMAX     = 5e-4        #  IDA maximum step[s]   ( hereAn adaptive timestepBDF DAE solver)

# Governor and stabilizer settings

const TG_TS  = 0.5    #Governor servo time constant[s]
const TG_EPS = 0.001   # other internal turbine stages are treated as very fast time (user can change)

const PSS_KS = 0.0    #PSS gain PSSFixed(V_pss = 0.0) made fixed 

#Droop settings
const R_SG_PU  = 0.05     #SG speed droop
const R_GFM_PU = 0.08    #GFM speed droop

#PLL settings
const GFL_PLL_KP = 0.02      #PLL Kp
const GFL_PLL_KI = 0.15        #PLL Ki

#Dynamic-line settings
const USE_SELECTIVE_DYNAMIC_LINES      = false     #false: this make all dynamic
const DYNAMIC_LINE_INCLUDE             = Tuple{Int,Int}[]    #Selected dynamic lines
const DYNAMIC_LINE_FORCE_STATIC        = Tuple{Int,Int}[]    #Lines kept static
const APPLY_DYNAMIC_LINE_DAMPING_EDITS = true     #Apply line R/X 
const DYNAMIC_LINE_R_SCALE             = 1.0    #Line resistance scale
const DYNAMIC_LINE_X_SCALE             = 1.0    #Line reactance scale
const DYNAMIC_LINE_R_FLOOR_PU          = 0.0005         #Minimum line resistance[p.u.]

#Small-signal report settings 
const RUN_SSA_REPORT                 = true     #Run small-signal report
const SSA_REPORT_LABEL               = "PSID_3BUS_EMT"   #Report label
const SSA_DAMPING_WARN_PCT           = 5.0    #Low damping threshold[%]
const SSA_REAL_NEAR_ZERO             = 1e-3   #Near-axis real-part limit
const SSA_FAST_REAL_THRESH           = 50.    #Fast-mode real-part threshold
const SSA_STABILITY_TOL              = 1e-8   #Stability tolerance
const SSA_TXT_FILENAME               = "ssa_stability_report.txt"#SSA text report
const SSA_EIG_CSV_FILENAME           = "ssa_eigenvalues.csv"#Eigenvalue CSV
const SSA_SUMMARY_CSV_FILENAME       = "ssa_summary_eigenvalues.csv"#Sorted eigenvalue CSV
const SSA_PLOT_FILENAME              = "ssa_eigs_real_vs_imag.png"#Eigenvalue plot
const SSA_FULL_MODAL_REPORT_FILENAME = "ssa_full_modal_report.txt"#Full modal report

# Resource partitions
const SG_BUSES      = [3]     #SG buses
const GFM_BUSES     = [1]  #GFM buses
const GFL_BUSES     = [2]#GFL buses
const SLACK_BUS_NUM = 3#Reference bus

#Export and plotting settings
const SAVE_DIAGNOSTIC_NETWORK_PQ = true   #Save network P/Q CSVs
const ABS_VOLTAGE_YLIM           = (0.94, 1.06)   #Voltage plot limits[p.u.]
const ABS_VOLTAGE_YTICKS         = 0.94:0.02:1.06  #Voltage plot ticks[p.u.]

const BUS_PLOT_COLORS_EXT = Dict{Int,String}(
    1 => "#1f77b4",  #Bus 1 colour
    2 => "#d62728",  #Bus 2 colour
    3 => "#2ca02c",   #Bus 3 colour
)

bus_plot_color(bus::Int) = get(BUS_PLOT_COLORS_EXT, bus, "#333333")

#Note to user: default parameters  maps are included  for all buses to make the script flexible when changing bus  assets assignments. 
#This prevents missing-parameter errors after changing bus assignments, and the defaults can be replaced with detailed values later.


#Per-bus inertia and damping
const H_DEFAULT    = 4.0     # Default inertia constant[s]
const D_SG_DEFAULT = 48.0   #Default SG damping coefficient

const H_ALL_map_ext = Dict{Int,Float64}(
    1 => 4.200,   #Bus 1 inertia constant[s]
    2 => 4.333,    #Bus 2 inertia constant[s]
    3 => 5.000,    #Bus 3 inertia constant[s]
)

const D_SG_map_ext = Dict{Int,Float64}(
    1 => D_SG_DEFAULT,     #Bus 1 SG damping coefficient
    2 => D_SG_DEFAULT,    #Bus 2 SG damping coefficient
    3 => 48.0      ,#Bus 3 SG damping coefficient
)

# GFM defaults and per-bus values
const GFM_TA_DEFAULT = 2.0 * 4.2      #  Virtual inertia time constant[s]
const GFM_KW_DEFAULT = 12.5#GFM frequency gain
const GFM_KD_DEFAULT = 12.5#GFM damping gain

const GFM_TA_map_ext = Dict{Int,Float64}(
    1 => 2.0 * 4.200,#Bus 1 GFM Ta[s]
    2 => 2.0 * 4.333,#Bus 2 GFM Ta[s]
    3 => 2.0 * 5.000,#Bus 3 GFM Ta[s]
)

const GFM_KW_map_ext = Dict{Int,Float64}(
    1 => GFM_KW_DEFAULT,#Bus 1 GFM frequency gain
    2 => GFM_KW_DEFAULT,#Bus 2 GFM frequency gain
    3 => GFM_KW_DEFAULT,#Bus 3 GFM frequency gain
)

const GFM_KD_map_ext = Dict{Int,Float64}(
    1 => GFM_KD_DEFAULT,#Bus 1 GFM damping gain
    2 => GFM_KD_DEFAULT,#Bus 2 GFM damping gain
    3 => GFM_KD_DEFAULT,#Bus 3 GFM damping gain
)



#Important Note to User  regarding GFL  kω_droop::Float64 = GFL_KW_DEFAULT
#I modified my local PSID source  to incoprate a FFR for my GFL outer power loop: eg 
# C:\Users\<my_username>\.julia\dev\PowerSimulationsDynamics\src\models\inverter_models\outer_control_models.jl
# The GFL ActivePowerPI/ReactivePowerPI mdl_outer_ode! block now includes FFR droop:
#   p_ref_eff = p_ref - Kω * (ω_pll - 1.0)
# and uses (p_ref_eff - p_oc) instead of the original fixed-P_ref error (p_ref - p_oc). User can go to the source file to make similkar chnages 
# Users keeping the original PSID source should not pass Kω/ext=Dict("Kω"=>...);
# use the default ActivePowerPI constructor 


# GFL defaults and per-bus values
const GFL_KW_DEFAULT = 10.0#GFL frequency-droop gain

const GFL_KW_map_ext = Dict{Int,Float64}(
    1 => GFL_KW_DEFAULT,#Bus 1 GFL frequency-droop gain
    2 => GFL_KW_DEFAULT,#Bus 2 GFL frequency-droop gain
    3 => GFL_KW_DEFAULT,#Bus 3 GFL frequency-droop gain
)

#Generator base-power map
const GEN_MBASE_DEFAULT = 300.0#Default generator base[MVA]

const GEN_MBASE_map_ext = Dict{Int,Float64}(
    1 => 300.0,#Bus 1 generator base[MVA]
    2 => 250.0,#Bus 2 generator base[MVA]
    3 => 400.0,#Bus 3 generator base[MVA]
)

getH(bus::Int)        = get(H_ALL_map_ext, bus, H_DEFAULT)
getDsg(bus::Int)      = get(D_SG_map_ext, bus, D_SG_DEFAULT)
getGfmTa(bus::Int)    = get(GFM_TA_map_ext, bus, GFM_TA_DEFAULT)
getGfmKw(bus::Int)    = get(GFM_KW_map_ext, bus, GFM_KW_DEFAULT)
getGfmKd(bus::Int)    = get(GFM_KD_map_ext, bus, GFM_KD_DEFAULT)
getGflKw(bus::Int)    = get(GFL_KW_map_ext, bus, GFL_KW_DEFAULT)
getGenMBase(bus::Int) = get(GEN_MBASE_map_ext, bus, GEN_MBASE_DEFAULT)