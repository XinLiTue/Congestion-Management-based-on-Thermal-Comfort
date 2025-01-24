module OptHP

using JuMP
using Interpolations
using Dates
using DataFrames
using Gurobi

export GEC, interpolate_data

# Constants used in the optimization problem
include("constants.jl")

# Utility functions (e.g. interpolation)
include("utils.jl")

# Main model 
include("grid.jl")
include("model.jl")

end
