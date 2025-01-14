module OptHP

using JuMP
using DataFrames
using Gurobi
using OffsetArrays

export GEC


include("heatpump.jl")
include("constants.jl")
include("user.jl")

end
