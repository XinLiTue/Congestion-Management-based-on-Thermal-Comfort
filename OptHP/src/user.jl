"""
Define users and their constraints.

Two types of users are defined:
- Users without a heat pump (RegularUser)
- Users with a heat pump (HeatPumpUser)

Both are subtypes of User.

The connectivity between users is defined by a tree. Some important characteristics
of trees are:
- Undirected graph in which any two vertices are connected by exactly one path
- No cycles, i.e. no closed loops
- Number of edges |E| is one less than the number of vertices |V|: |E| = |V| - 1

Each user has the following characteristics:
- Identified by unique user ID (uid)
- Has exactly 1 parent node (except the root node, also called the slack bus)
- Has a set of child nodes
"""


abstract type Bus end

# the slack bus is the root node of the tree, it has no parent
struct SlackBus <: Bus
    node::Int
    children::Vector{<:Bus}
end

# regular user bus without a heat pump
struct RegularBus <: Bus
    node::Int
    parent::Int
    children::Vector{Int}
end

# user bus with a heat pump
struct HeatPumpBus <: Bus
    node::Int
    parent::Int
    children::Vector{Int}
end

# a single line connects two users, the parent bus and the child bus
struct Line
    parent::Int
    child::Int
    length::Float64
    cable::String
    Inom::Float64
    R::Float64
    X::Float64
end

function Base.show(io::IO, b::Bus)
    println(io, "Bus: $(b.node) with parent: $(b.parent) and children: $(b.children)")
end

function Base.show(io::IO, l::Line)
    println(io, "Line: $(l.parent) -> $(l.child) [X: $(l.X), R: $(l.R), Inom: $(l.Inom)]")
end