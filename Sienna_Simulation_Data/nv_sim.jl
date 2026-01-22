using Pkg
script_path = @__DIR__
Pkg.activate(script_path)
using Logging
using PowerSystems
const PSY = PowerSystems
using PowerSimulations
using Xpress
using Gurobi
const PSI = PowerSimulations
using InfrastructureSystems
const IS = InfrastructureSystems
using HydroPowerSimulations
using StorageSystemsSimulations
const SSS = StorageSystemsSimulations
const HPS = HydroPowerSimulations
using PowerNetworkMatrices
const PNM = PowerNetworkMatrices
using Dates
using TimeSeries
using Random
using CSV
using DataFrames
using Dates
using Plots
using DataStructures
using JSON3
using DelimitedFiles
using JLD2
using Logging
using HDF5
using HiGHS

sys = System(joinpath(script_path, "..", "NV_sys_new.json"))

# Modifying reactance (line.x) and resistance (line.r) of the lines since it cannot be zero 

map( line -> (if (line.r == 0.0 || line.x == 0.0)
        @info("Fixing Line: ", line.name)
        set_r!(line, Float64(0.1))  # Assign default value for r = resistance
        set_x!(line, Float64(0.1))  # Assign default value for x = reactance
    end
            ),  get_components(Line, sys))


mip_gap = 0.0125 

#=
optimizer = optimizer_with_attributes(
    Xpress.Optimizer,
    "MIPRELSTOP" => 1e-4,   # relative MIP gap
)=#

optimizer = optimizer_with_attributes(
                HiGHS.Optimizer,
                #"parallel" => "on",
                "mip_rel_gap" => mip_gap)

transform_single_time_series!(sys, Hour(24), Day(1))

template_uc =
    ProblemTemplate(
        NetworkModel(CopperPlatePowerModel;
        use_slacks = true,
        ),
    )
#=
 ## Remove/Disable interchange areas    
for ai in get_components(AreaInterchange, sys) 
 
    #set_available!(ai, false)
    remove_component!(sys, ai)
end
=#

     

set_device_model!(template_uc, ThermalStandard, ThermalStandardUnitCommitment)
set_device_model!(template_uc, ThermalMultiStart, ThermalStandardUnitCommitment)
#set_device_model!(template_uc, ThermalStandard, ThermalStandardUnitCommitment)
#set_device_model!(template_uc, ThermalMultiStart, ThermalStandardUnitCommitment)
set_device_model!(template_uc, RenewableDispatch, RenewableFullDispatch)
set_device_model!(template_uc, RenewableNonDispatch, FixedOutput)
set_device_model!(template_uc, PowerLoad, StaticPowerLoad)
set_device_model!(template_uc, HydroDispatch, HydroDispatchRunOfRiver) #HydroDispatchRunOfRiver


model = DecisionModel(
    template_uc,
    sys;
    name = "UC",
    optimizer = optimizer,
    system_to_file = false,
    initialize_model = true,
    check_numerical_bounds = false,
    optimizer_solve_log_print = true,
    direct_mode_optimizer = false,
    rebuild_model = false,
    store_variable_names = true,
    calculate_conflict = true,
)

models = SimulationModels(
            decision_models = [model],
        )

DA_sequence = SimulationSequence(
    models = models,
    ini_cond_chronology = InterProblemChronology(),
)    

initial_date = "2019-01-01"
steps_sim    = 3
current_date = string( today() )
sim = Simulation(
    name = current_date * "_NV" * "_" * string(steps_sim) * " steps",
    steps = steps_sim,
    models = models,
    initial_time = DateTime("2030-01-01T00:00:00"),
    sequence = DA_sequence,
    simulation_folder = "."#".",  use "tempdir()" if you dont want to store simulation data
)

build!(sim) 

execute!(sim)