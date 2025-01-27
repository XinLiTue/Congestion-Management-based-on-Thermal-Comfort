module OptHP

using JuMP
using Interpolations
using Dates
using DataFrames
using Gurobi
using Clarabel

export GEC, interpolate_data
export S_base, V_base, Z_base, I_base

# Constants used in the optimization problem
include("constants.jl")

# Utility functions (e.g. interpolation)
include("utils.jl")

# Main model 
include("grid.jl")
include("model.jl")

end
