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
using PlotlyJS
using PowerGraphics
using YAML

sys = System(joinpath(script_path, "..", "NV_sys_new.json"))
# Update the hydro
# set limit and set active power limits
# Check the capacity for HydroDispatch. If 120 MW, change to 1.2 GWh so 10X.

hydro = collect(get_components(HydroDispatch, sys))[1]
hydro.active_power_limits = (min = 0.0, max = 1210.0)
@info "Current max power: ", hydro.active_power_limits.max, " MW"
@info "Current time limit: ", hydro.time_limits
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
steps_sim    = 365
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

# Reading results
results = SimulationResults(sim)

uc_results = get_decision_problem_results(results, "UC"); # UC stage result metadata

read_variables(uc_results) # read all the result variables
read_parameters(uc_results) # read all the parameters
list_variable_names(uc_results) #  list the variable names contained in uc_results
list_parameter_names(uc_results) # list the parameters contained in uc_results

# Plotting #TODO
#=
# 1. Graph that shows the total installed capacity of the system and broken down by what kind of generation.  Hydro, solar, NG, biomass, petroleum
Simulation is not needed for this.
# 2. Graph that shows the total GWh of the system. All energy produced during the year and breaking down by type of generation. TOTALS
Refer to Power Graphics. 
Function plotfuel() sort of. Generator mapping file is needed.
3. Plot of the load across the year. 
read_realized_parameter(uc_results, "ActivePowerTimeSeriesParameter__PowerLoad", table_format = TableFormat.WIDE)
Multiply by -1 to make it positive.
4. Make a pptx for Anna to add the hydro
5. Check the capacity of the system 1.2 GW but now is 120 MW and if so, we need to add 
=#

# 1. Graph that shows the total installed capacity of the system and broken down by what kind of generation.  Hydro, solar, NG, biomass, petroleum

palette = YAML.load_file("color-palette.yaml")

thermals = collect(get_components(ThermalStandard, sys))
hydros   = collect(get_components(HydroDispatch, sys))
renewable_dispatch = collect(get_components(RenewableDispatch, sys))
renewable_nondispatch = collect(get_components(RenewableNonDispatch, sys))

unique([g.fuel for g in thermals])
unique(g.prime_mover_type for g in renewable_dispatch)
unique(g.prime_mover_type for g in renewable_nondispatch)


capacity = Dict{String, Float64}()

capacity["Hydro"] = sum(g.active_power_limits.max for g in hydros)
capacity["Solar_RD"] = sum( g.rating for g in renewable_dispatch if g.prime_mover_type == PrimeMovers.PVe) # Is rating correct? there is no active_power_limits
capacity["Solar_ND"] = sum( g.rating for g in renewable_nondispatch if g.prime_mover_type == PrimeMovers.PVe)  # Is rating correct? there is no active_power_limits
capacity["Wind_RD"] = sum( g.rating for g in renewable_dispatch if g.prime_mover_type == PrimeMovers.WT)  # Is rating correct? there is no active_power_limits
capacity["Natural Gas"] = sum( g.active_power_limits.max for g in thermals if g.fuel == ThermalFuels.NATURAL_GAS)
capacity["Geothermal"] = sum(g.rating for g in thermals if g.fuel == ThermalFuels.GEOTHERMAL)
capacity["Distillate_Fuel_Oil"] = sum(g.rating for g in thermals if g.fuel == ThermalFuels.DISTILLATE_FUEL_OIL)
capacity["Wood"] = sum(g.rating for g in thermals if g.fuel == ThermalFuels.WOOD_WASTE_SOLIDS)


labels = collect(keys(capacity))

mw = [capacity[label] for label in labels]

# 8 distinct colors
colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#8c564b", 
          "#d62728", "#9467bd", "#17becf", "#e377c2"]

trace = PlotlyJS.bar(
x = labels,
y = mw,
marker_color = colors
)

layout = PlotlyJS.Layout(
title = "Installed Generation Capacity by Technology",
yaxis_title = "Installed Capacity (MW)",
xaxis_title = "Generation Type"
)

PlotlyJS.plot(trace, layout)

#### Another way with PowerGraphics ##############################

using PlotlyJS, YAML


# Load palette YAML

palette_file = "color-palette.yaml"
palette_data = YAML.load_file(palette_file)

# Convert to Dict{String,String} for easy lookup
palette_dict = Dict{String,String}()
for (k,v) in palette_data
    palette_dict[k] = v["RGB"]
end


# Collect components

thermals = collect(get_components(ThermalStandard, sys))
hydros   = collect(get_components(HydroDispatch, sys))
renewable_dispatch = collect(get_components(RenewableDispatch, sys))
renewable_nondispatch = collect(get_components(RenewableNonDispatch, sys))


# Compute capacities with correct palette keys

capacity = Dict{String,Float64}()
colors   = Dict{String,String}()

# Thermals
for g in thermals
    key = nothing
    if g.fuel == ThermalFuels.NATURAL_GAS && g.prime_mover_type == PrimeMovers.CT
        key = "Natural Gas CT"
    elseif g.fuel == ThermalFuels.NATURAL_GAS && g.prime_mover_type == PrimeMovers.ST
        key = "Natural Gas ST"
    elseif g.fuel == ThermalFuels.NATURAL_GAS && g.prime_mover_type == PrimeMovers.CC
        key = "Natural Gas CC"
    elseif g.fuel == ThermalFuels.DISTILLATE_FUEL_OIL
        key = "Petroleum"
    elseif g.fuel == ThermalFuels.GEOTHERMAL
        key = "Geothermal"
    elseif g.fuel == ThermalFuels.WOOD_WASTE_SOLIDS
        key = "Biopower"
    else
        key = string(g.fuel)
    end

    capacity[key] = get(capacity, key, 0.0) + g.active_power_limits.max
    colors[key] = get(palette_dict, key, "grey")
end

# Hydros
for g in hydros
    key = "Hydropower"
    capacity[key] = get(capacity, key, 0.0) + g.active_power_limits.max
    colors[key] = get(palette_dict, key, "grey")
end

# Renewables Dispatchable
for g in renewable_dispatch
    key = nothing
    if g.prime_mover_type == PrimeMovers.PVe
        key = "PV"
    elseif g.prime_mover_type == PrimeMovers.WT
        key = "Wind"
    else
        key = string(g.prime_mover_type)
    end
    capacity[key] = get(capacity, key, 0.0) + g.rating
    colors[key] = get(palette_dict, key, "grey")
end

# Renewables Non-Dispatchable
for g in renewable_nondispatch
    key = nothing
    if g.prime_mover_type == PrimeMovers.PVe
        key = "PV"
    else
        key = string(g.prime_mover_type)
    end
    capacity[key] = get(capacity, key, 0.0) + g.rating
    colors[key] = get(palette_dict, key, "grey")
end


# PlotlyJS bar

labels = collect(keys(capacity))
mw     = [capacity[l] for l in labels]
color_list = [colors[l] for l in labels]

trace = PlotlyJS.bar(
    x = labels,
    y = mw,
    marker_color = color_list
)

layout = PlotlyJS.Layout(
    title = "Installed Generation Capacity by Technology",
    yaxis_title = "Installed Capacity (MW)",
    xaxis_title = "Generation Type"
)

PlotlyJS.plot(trace, layout)


# 2. Graph that shows the total GWh of the system. All energy produced during the year and breaking down by type of generation.

palette = load_palette("color-palette.yaml")

plot = plot_fuel(
    uc_results;
    generator_mapping_file = "generator_mapping.yaml",
    bar = true,
    stack = true,
    title = "Total Energy Produced in 2030 by Fuel",
    palette = palette, 
    set_display = true,
    variables = [:generation],
    slacks = false,
    curtailment = false,
    load = false,
    y_label = "Energy (MWh)"
)

#=
3. Plot of the load across the year. 
read_realized_parameter(uc_results, "ActivePowerTimeSeriesParameter__PowerLoad", table_format = TableFormat.WIDE)
Multiply by -1 to make it positive.
=#


# Read load data and make it positive

load_df = read_realized_parameter(
    uc_results, 
    "ActivePowerTimeSeriesParameter__PowerLoad", 
    table_format = TableFormat.WIDE
)

for col in names(load_df)[2:end]  # skip DateTime
    load_df[!, col] .= -1 .* load_df[!, col]
end




# Plot:

traces = PlotlyJS.GenericTrace[]  # GenericTrace

for col in names(load_df)[2:end]
    push!(traces, PlotlyJS.scatter(
        x = load_df.DateTime,
        y = load_df[!, col],
        mode = "lines",
        name = col
    ))
end


layout = PlotlyJS.Layout(
    title = "Load Across the Year",
    xaxis_title = "Time",
    yaxis_title = "Load (MW)",
    hovermode = "x unified"
)

plt = PlotlyJS.Plot(traces, layout)
display(plt)
