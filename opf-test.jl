using JuMP
using Gurobi
using Clarabel
using CairoMakie
using Random


Random.seed!(333)


function plot_line_currents(model::Model, lines::Vector{Tuple{Int,Int}}=[(34, 0), (2, 74)])
    fig = Figure(; size=(1000, 600))

    I_line = value.(model[:I_line])
    for (idx, (i, j)) in enumerate(lines)
        ax = Axis(fig[idx, 1], xlabel="Time [hours]",
            ylabel="Current [A]",
            title="($i -> $j)",
            xticks=(1:4:97, string.(0:1:24))
        )
        scatterlines!(ax, Vector(I_line[(i, j), :]), color=:blue, label="$i->$j", linewidth=2)
    end

    fig
end

function plot_pq(model::Model, bus::Int)
    fig = Figure(; size=(1000, 300))
    ax = Axis(fig[1, 1], xlabel="Time [hours]",
        ylabel="Power [p.u.]",
        title="Power flow",
        xticks=(1:4:97, string.(0:1:24))
    )

    # index of slack bus (transformer)
    # SB = argmin(value.(model[:P]).axes[2])
    P_trafo = Matrix(value.(model[:P]))[bus, :]
    Q_trafo = Matrix(value.(model[:Q]))[bus, :]


    # lines!(ax, sol[:P], color = :blue, label = "Transformer", linestyle = :so
    scatterlines!(ax, P_trafo, color=:blue, label="P trafo", linewidth=2)
    scatterlines!(ax, Q_trafo, color=:red, label="Q trafo", linewidth=2)
    axislegend(ax, position=:lb)

    # display
    fig
end



# p.u. base values
S_base = 1E5 # [VA]
V_base = 230.0 # [V]

Z_base = V_base^2 / S_base
I_base = S_base / V_base
S_max = 1E6 * 0.23 * 0.16 * 3 / S_base

# voltage constraints V_base ± 10% -> (207V, 253V) at 230V base
V_lb = 0.90 # [pu]
V_ub = 1.10

T = 1:96


function add_grid_model(;
    model::Model,
    T=1:96
)
    # test network with 4 buses 
    #       0   (slack)
    #       |
    #       1   (junction)
    #      / \
    #     2   3 (loads)

    B = [0, 1, 2, 3]
    H = [2, 3]
    L = [(0, 1), (1, 2), (1, 3)]

    bus_out = Dict(0 => [1], 1 => [2, 3], 2 => [], 3 => [])
    bus_in = Dict(0 => [], 1 => [0], 2 => [1], 3 => [1])

    # base loads 
    # P_base = 1e4 * ones(length(T)) / S_base
    # P_base = round.(1e3 .* abs.(randn(length(T))) / S_base .+ 3e4 / S_base, digits=3)
    P_base = 1e3 .* abs.(randn(length(T))) / S_base .+ 2e4 / S_base
    println(P_base)
    Q_base = zeros(length(T))

    # slack bus ID
    SB = 0

    # real and reactive line impedances
    R = Dict((i, j) => 0.1 / Z_base for (i, j) in L)
    X = Dict((i, j) => 0.1 / Z_base for (i, j) in L)
    Inom = Dict((i, j) => 1.0 for (i, j) in L)

    # variables
    @variables(model, begin
        # voltage squared eq. (7)
        V_lb^2 <= v[B, T] <= V_ub^2

        # bus power injections
        P[B, T], (base_name = "PBusInjection")
        Q[B, T], (base_name = "QBusInjection")

        # power flows
        P_line[L, T], (base_name = "PLineFlow")
        Q_line[L, T], (base_name = "QLineFlow")

        # line current
        0 <= I_line[L, T], (base_name = "CurrentSquare")
    end)

    # slack bus constraints
    @constraints(model, begin
        TrafoLimit[t in T], [S_max, P[SB, t], Q[SB, t]] in SecondOrderCone()
    end)

    # bus / line constraints
    @constraints(model, begin
        # real power balance eq. (1)
        RealPowerBalance[j in B, t in T],
        P[j, t] == sum(P_line[(j, k), t] for k in bus_out[j]) -
                   sum(P_line[(i, j), t] - R[(i, j)] * I_line[(i, j), t] for i in bus_in[j])

        # reactive power balance eq. (2)
        ReactivePowerBalance[j in B, t in T],
        Q[j, t] == sum(Q_line[(j, k), t] for k in bus_out[j]) -
                   sum(Q_line[(i, j), t] - X[(i, j)] * I_line[(i, j), t] for i in bus_in[j])

        # voltage drop eq. (3)
        VoltageDrop[(i, j) in L, t in T],
        v[j, t] == v[i, t] - 2 * (P_line[(i, j), t] * R[(i, j)] + Q_line[(i, j), t] * X[(i, j)]) +
                   (R[(i, j)]^2 + X[(i, j)]^2) * I_line[(i, j), t]

        # conic OPF eq. (4)
        ConicOPF[(i, j) in L, t in T],
        [
            I_line[(i, j), t] + v[i, t],
            2 * P_line[(i, j), t],
            2 * Q_line[(i, j), t],
            I_line[(i, j), t] - v[i, t],
        ] in SecondOrderCone()

        # line current limit eq. (6)
        LineCurrentLimit[(i, j) in L, t in T], I_line[(i, j), t] <= Inom[(i, j)]
    end)

    # load constraints
    @constraints(model, begin
        Transmission[i in B, t in T; i ∉ H && i ≠ SB], P[i, t] == 0
        RealBaseLoad[i in H, t in T], P[i, t] == -P_base[t]
        ReactiveBaseLoad[i in H, t in T], Q[i, t] == -Q_base[t]
    end)
end


# create the model
# model = Model(Clarabel.Optimizer)
model = Model(Gurobi.Optimizer)
set_attribute(model, "BarHomogeneous", 1)
set_attribute(model, "Presolve", 0)

# add grid model 
add_grid_model(model=model, T=T)

P = model[:P]

# define objective functions
@expressions(model, begin
    J_gen, sum(P[0, T])
end)

@objective(model, Min, J_gen)

# optimize the model and return it
optimize!(model)
@assert is_solved_and_feasible(model)
