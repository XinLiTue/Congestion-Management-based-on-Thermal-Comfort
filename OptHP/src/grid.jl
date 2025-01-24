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
    PV::Float64         # PV capacity [MW]
    A::Matrix{Float64}
    B::Matrix{Float64}
    T_supply::Float64   # supply temperature [°C]

    function HeatPumpBus(
        node::Int,
        adjacent::Set{Int},
        PV::Float64,
        A::Matrix{Float64},
        B::Matrix{Float64},
        T_supply::Float64=40.0
    )
        if size(A) != (3, 3) || size(B) != (3, 3)
            throw(ArgumentError("A and B must be 3x3 matrices"))
        end
        new(node, adjacent, PV, A, B, T_supply)
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
    T::UnitRange{Int}                       # Set of discrete time steps
    Δt::Float64                             # Time step duration

    # constructor with default time step duration
    function Grid(
        buses::NTuple{N,AbstractBus},
        lines::NTuple{M,Line},
        T::UnitRange{Int},
        Δt::Float64=1.0) where {N,M}

        new(buses, lines, T, Δt)
    end
end

# build the grid network consisting of buses and lines, 
# defined by the network and connections dataframes
function Grid(
    network::DataFrame,
    connections::DataFrame,
    meta::Dict,
    T::UnitRange{Int}
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
        PV = row.PV * 1e-3

        if row.HP == 1
            users[node] = Bus(node, adjacent, PV, meta["H14"]["A"], meta["H14"]["B"])
        else
            users[node] = Bus(node, adjacent, PV)
        end
    end

    # create grid
    lines = [Line(users[row.StartNode], users[row.EndNode],
        row.Length, 1E3 * row.Inom, row.R, row.X) for row in eachrow(network)]

    grid = Grid(Tuple(values(users)), Tuple(lines), T)
    return grid
end

# get the sets of buses and users
function get_sets(grid::Grid)
    B = Set(getfield.(grid.buses, :node))
    H = Set([bus.node for bus in grid.buses if bus isa User])
    H_HP = Set([bus.node for bus in grid.buses if bus isa HeatPumpBus])
    notH = setdiff(B, H)
    return (B=B, H=H, H_HP=H_HP, notH=notH, T=grid.T)
end

# get the set of lines indexed by (start, stop) bus pairs
function get_line_set(grid::Grid)
    return Set([(line.start.node, line.stop.node) for line in grid.lines])
end

# get the buses going into or outgoing of this bus
function bus_in(bus::Int, grid::Grid)
    return Set([line.start.node for line in grid.lines if line.stop.node == bus])
end
function bus_out(bus::Int, grid::Grid)
    return Set([line.stop.node for line in grid.lines if line.start.node == bus])
end

# get lines going into or out of this bus
function get_incoming_lines(bus::Int, grid::Grid)
    return [line for line in grid.lines if line.stop.node == bus]
end
function get_outgoing_lines(bus::Int, grid::Grid)
    return [line for line in grid.lines if line.start.node == bus]
end

# show methods for pretty printing
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
    println(io, "Grid with $(length(g.buses)) buses and $(length(g.lines)) lines over $(length(g.T)) time steps")
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

# get the slack bus and verify it's unique
function get_slack_bus(grid::Grid)
    slack_buses = [bus for bus in grid.buses if bus isa SlackBus]
    @assert length(slack_buses) == 1
    return slack_buses[1]
end

# get the user buses without heat pumps
function get_nonhp_buses(grid::Grid)
    return [bus for bus in grid.buses if bus isa UserBus]
end

# get the user buses with heat pumps
function get_hp_buses(grid::Grid)
    return [bus for bus in grid.buses if bus isa HeatPumpBus]
end

# get user buses
function get_user_buses(grid::Grid)
    return [bus for bus in grid.buses if bus isa User]
end