module OptHP

using JuMP
using DataFrames
using Gurobi

export GEC


include("heatpump.jl")
include("constants.jl")
include("grid.jl")
include("model.jl")

end
