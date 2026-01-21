data_dir = @__DIR__
using CSV 
using DataFrames
### Code to Update Hydro Modeling 
unit_dict = Dict(
    "N7" => 127.0,
    "N1" => 130.0, 
    "N2" => 130.0, 
    "N3" => 130.0, 
    "N4" => 130.0,
    "N5" => 130.0,
    "N6" => 130.0,
    "N8" => 130.0,
)

# ---------------------
#    Reservoir Data
# ---------------------

min_res_level = 950 #ft
dam_volume = 118.8e6 #cubic ft
res_levels = CSV.read(joinpath(data_dir, "historic_res_levels.csv"), DataFrame)

