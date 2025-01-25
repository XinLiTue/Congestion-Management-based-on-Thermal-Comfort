"""Define model related functions."""

function add_grid_model(;
    model::Model,
    grid::Grid,
    limit::Tuple,
    loads_real::DataFrame,
    loads_reactive::DataFrame
)
    # sets
    sets = get_sets(grid)
    B, H, H_HP, T = sets.B, sets.H, sets.H_HP, sets.T
    L = get_line_set(grid)

    I_max = 1e4 # [A]


    # base loads 
    P_base = loads_real
    Q_base = loads_reactive

    # slack bus ID
    SB = get_slack_bus(grid).node

    # real and reactive line impedances
    R = Dict((line.start.node, line.stop.node) => line.R for line in grid.lines)
    X = Dict((line.start.node, line.stop.node) => line.X for line in grid.lines)

    # variables
    @variables(model, begin
        # voltage squared eq. (7)
        (V_lb * V_ref)^2 <= v[B, T] <= (V_ub * V_ref)^2

        # bus power injections
        P[B, T], (base_name = "PBusInjection")
        Q[B, T], (base_name = "QBusInjection")

        # power flows
        P_line[L, T], (base_name = "PLineFlow")
        Q_line[L, T], (base_name = "QLineFlow")

        # line current
        0 <= I_line[L, T], (base_name = "CurrentSquare")

        # heat pump electrical power consumption
        # P_HP[H_HP, T], (base_name = "HPElectricalPower")
    end)

    # slack bus constraints
    @constraints(model, begin
        TrafoPowerLimitForCongestion[t in limit[1]], P[SB, t] <= limit[2]
        TrafoLimit[t in T], [S_max, P[SB, t], Q[SB, t]] in SecondOrderCone()
    end)

    # bus / line constraints
    @constraints(model, begin
        # real power balance eq. (1)
        RealPowerBalance[j in B, t in T],
        P[j, t] == sum(P_line[(j, k), t] for k in bus_out(j, grid)) -
                   sum(P_line[(i, j), t] - R[(i, j)] * I_line[(i, j), t] for i in bus_in(j, grid))

        # reactive power balance eq. (2)
        ReactivePowerBalance[j in B, t in T],
        Q[j, t] == sum(Q_line[(j, k), t] for k in bus_out(j, grid)) -
                   sum(Q_line[(i, j), t] - X[(i, j)] * I_line[(i, j), t] for i in bus_in(j, grid))

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
        # ConicOPF[(i, j) in L, t in T],
        # I_line[(i, j), t] >= (P_line[(i, j), t]^2 + Q_line[(i, j), t]^2) / v[i, t]


        # line current limit eq. (6)
        LineCurrentLimit[(i, j) in L, t in T], I_line[(i, j), t] <= I_max
    end)

    # load constraints
    @constraints(model, begin
        Transmission[i in B, t in T; i ∉ H && i ≠ SB], P[i, t] == 0
        RealBaseLoad[i in H, t in T], P[i, t] == -P_base[t, "$i"]
        # [i in H, t in T; i ∉ H_HP], P[i, t] == -P_base[t, "$i"]
        # [i in H_HP, t in T], P[i, t] == -P_base[t, "$i"] - P_HP[i, t]
        ReactiveBaseLoad[i in H, t in T], Q[i, t] == -Q_base[t, "$i"]
    end)

    # ### TESTING ###
    # @constraints(model, begin
    #     # for testing we set P_hp to zero
    #     [i in H_HP, t in T], P_HP[i, t] == 0
    # end)
    # ### TESTING ###
end


function GEC(;
    network::DataFrame,
    connections::DataFrame,
    loads_real::DataFrame,
    loads_reactive::DataFrame,
    weather::DataFrame,
    limit::Tuple=(1:96, 250 * 1e3),
    meta::Dict=Dict(),
    silent=true,
    T=1:96
)
    # create the model
    model = Model(Gurobi.Optimizer)
    if silent
        set_silent(model)
    end

    # apply T to all inputs
    loads_real = loads_real[T, :]
    loads_reactive = loads_reactive[T, :]
    weather = weather[T, :]


    # add grid model 
    grid = Grid(network, connections, meta, T)
    add_grid_model(model=model,
        grid=grid,
        limit=limit,
        loads_real=loads_real,
        loads_reactive=loads_reactive,
    )
    # println("grid.lines: ", grid.lines)

    # variables
    I_line = model[:I_line]
    P = model[:P]

    # define objective functions
    @expressions(model, begin
        J_loss, sum(l.R * I_line[(l.start.node, l.stop.node), t] for t in grid.T, l in grid.lines)
        J_gen, sum(P[0, T]) * 1e3
    end)
    @objective(model, Min, J_gen)
    # set_attribute(model, "BarHomogeneous", 1)
    # set_attribute(model, "NumericFocus", 3)
    # MOI.set(model, MOI.RelativeGapTolerance(), 1E-3)
    set_optimizer_attribute(model, "LogFile", "my_log_file.txt")


    # optimize the model and return it
    optimize!(model)
    # @assert is_solved_and_feasible(model)
    return model
end
