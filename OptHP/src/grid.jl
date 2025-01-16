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


abstract type AbstractBus end
abstract type User <: AbstractBus end

# the slack bus is the root node of the tree, it has no parent
struct SlackBus <: AbstractBus
    node::Int
    adjacent::Set{Int}
end

# junction bus without a user, connects other buses
struct JunctionBus <: AbstractBus
    node::Int
    adjacent::Set{Int}

    function JunctionBus(node::Int, adjacent::Set{Int})
        if node == 0
            return SlackBus(node, adjacent)
        end
        new(node, adjacent)
    end
end

# regular user bus without a heat pump
struct UserBus <: User
    node::Int
    adjacent::Set{Int}
    PV::Float64
end

# user bus with a heat pump
struct HeatPumpBus <: User
    node::Int
    adjacent::Set{Int}
    PV::Float64
    A::Matrix{Float64}
    B::Matrix{Float64}

    function HeatPumpBus(node::Int, adjacent::Set{Int}, PV::Float64, A::Matrix{Float64}, B::Matrix{Float64})
        if size(A) != (3, 3) || size(B) != (3, 3)
            throw(ArgumentError("A and B must be 3x3 matrices"))
        end
        new(node, adjacent, PV, A, B)
    end
end

# outer constructors for Bus structs
Bus(node::Int, adjacent::Set{Int}) = JunctionBus(node, adjacent)
Bus(node::Int, adjacent::Vector{Int}) = JunctionBus(node, Set(adjacent))
Bus(node::Int, adjacent::Set{Int}, PV::Float64) = UserBus(node, adjacent, PV)
Bus(node::Int, adjacent::Vector{Int}, PV::Float64) = UserBus(node, Set(adjacent), PV)
Bus(node::Int, adjacent::Set{Int}, PV::Float64, A::Matrix{Float64}, B::Matrix{Float64}) =
    HeatPumpBus(node, adjacent, PV, A, B)

# a single line connects two buses
struct Line
    start::AbstractBus
    stop::AbstractBus
    length::Float64
    Inom::Float64
    R::Float64
    X::Float64
end

# a grid network struct 
struct Grid
    buses::NTuple{N,AbstractBus} where {N}  # Tuple with N elements of AbstractBus
    lines::NTuple{M,Line} where {M}         # Tuple with M elements of Line
end


function Base.show(io::IO, b::AbstractBus)
    println(io, "$(split(string(typeof(b)),'.')[2]): $(b.node) with adjacent: [$(join(b.adjacent|>collect|>sort,','))]")
end

function Base.show(io::IO, b::UserBus)
    println(io, "UserBus: $(b.node) with PV: $(b.PV) and adjacent: [$(join(b.adjacent|>collect|>sort,','))]")
end

function Base.show(io::IO, l::Line)
    println(io, "$(l.start.node) -> $(l.stop.node) [X: $(l.X), R: $(l.R), Inom: $(l.Inom)]")
end

function Base.show(io::IO, g::Grid)
    println(io, "Grid with $(length(g.buses)) buses and $(length(g.lines)) lines")
end

# figure out which buses are connected to this bus
function connected_buses(bus::Int, lines::Set{Tuple{Int64,Int64}})
    connected = Set{Int}()
    for line in lines
        for idx in [1, 2]
            if line[idx] == bus
                push!(connected, line[3-idx])
            end
        end
    end
    return connected
end



function build_grid(
    network::DataFrame,
    connections::DataFrame,
    meta::Dict
)
    # a graph G = (V, E) is made of buses (V) and lines (E)
    V = Set([network.StartNode; network.EndNode])
    E = Set([(row.StartNode, row.EndNode) for row in eachrow(network)])
    @assert length(E) == length(V) - 1

    # sets
    B = Set(0:nrow(network))                                # all buses [-]
    H = Set(connections[!, :Node])                          # user buses [-]
    notH = setdiff(B, H)                                    # non-user buses [-]
    H_HP = Set(connections[connections.HP.==1, :Node])      # buses with HP [-]
    @assert issubset(H_HP, H)

    users = Dict{Int,AbstractBus}()

    # non-user bus connections
    for node in notH
        adjacent = connected_buses(node, E)
        users[node] = Bus(node, adjacent)
    end

    # user bus connections
    for row in eachrow(connections)
        node = row.Node
        adjacent = connected_buses(node, E)

        if row.HP == 1
            users[node] = Bus(node, adjacent, row.PV, meta["H14"]["A"], meta["H14"]["B"])
        else
            users[node] = Bus(node, adjacent, row.PV)
        end
    end

    # create grid
    lines = [Line(users[row.StartNode], users[row.EndNode],
        row.Length, row.Inom, row.R, row.X) for row in eachrow(network)]

    grid = Grid(Tuple(values(users)), Tuple(lines))
    return grid
end