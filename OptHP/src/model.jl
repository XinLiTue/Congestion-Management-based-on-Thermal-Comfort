"""Define model related functions."""

function add_grid_model(;
    model::Model,
    grid::Grid,
    limit::Tuple,
    df::DataFrame,
    medians::DataFrame,
    model_hp::Bool=false
)
    # sets
    sets = get_sets(grid)
    B, H, H_HP, T = sets.B, sets.H, sets.H_HP, sets.T
    L = get_line_set(grid)

    # base loads 
    P_base = df.usage_total_pu
    Q_base = zeros(T)

    # slack bus ID
    SB = get_slack_bus(grid).node

    # real and reactive line impedances
    R = Dict((line.start.node, line.stop.node) => line.R for line in grid.lines)
    X = Dict((line.start.node, line.stop.node) => line.X for line in grid.lines)
    Inom = Dict((line.start.node, line.stop.node) => line.Inom for line in grid.lines)

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

    # # slack bus constraints
    # @constraints(model, begin
    #     TrafoPowerLimitForCongestion[t in limit[1]], P[SB, t] <= limit[2]
    #     # TrafoLimit[t in T], [S_max, P[SB, t], Q[SB, t]] in SecondOrderCone()
    # end)

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

        # line current limit eq. (6)
        # LineCurrentLimit[(i, j) in L, t in T], I_line[(i, j), t] <= Inom[(i, j)]
    end)

    # load constraints
    @constraints(model, begin
        Transmission[i in B, t in T; i ∉ H && i ≠ SB], P[i, t] == 0

        # base load constraints
        RealBaseLoad[i in H, t in T; i ∉ H_HP], P[i, t] == -P_base[t]
        ReactiveBaseLoad[i in H, t in T; i ∉ H_HP], Q[i, t] == -Q_base[t]

        # grid congestion constraint 
        # 1 heat pump consumes 1.8 kW electric
        # in p.u. this is 1.8 / S_base = 1.8 / 1E5 = 1.8E-5
        # CongestionLimit[t in limit[1]], sum(P_HP[:, t]) <= limit[2]
    end)

    if model_hp
        @variables(model, begin
            P_HP[H_HP, T], (base_name = "HPElectricalPower")
            Q_HP[H_HP, T], (base_name = "HPReactivePower")
        end)

        @constraints(model, begin
            RealHPBaseLoad[i in H_HP, t in T], P[i, t] == -P_base[t] - P_HP[i, t]
            ReactiveHPBaseLoad[i in H_HP, t in T], Q[i, t] == -Q_base[t] - Q_HP[i, t]
            QvsP[i in H_HP, t in T], Q_HP[i, t] == tan_phi_load * P_HP[i, t]
        end)

        # add heat pump model
        add_heatpump_model(model, grid, df)
        add_ppd(model, grid, medians)

        @expressions(model, begin
            P_hp_user_load[t in T], sum(P_HP[i, t] for i in H_HP)
            P_non_hp_user_load[t in T], sum(P[i, t] for i in H if i ∉ H_HP)
            J_heat, sum(J_heat) # in €
            J_ppd, sum(0.27 * PPD[:, :]) # in %
        end)
    else
        @constraints(model, begin
            RealHPBaseLoad[i in H_HP, t in T], P[i, t] == -P_base[t]
            ReactiveHPBaseLoad[i in H_HP, t in T], Q[i, t] == -Q_base[t]
        end)
    end
end


function GEC(;
    network::DataFrame,
    connections::DataFrame,
    df::DataFrame,
    medians::DataFrame,
    limit::Tuple=(40:48, 0.0),
    meta::Dict=Dict(),
    silent=true,
    T=1:96,
    model_hp::Bool=false
)
    # create the model
    model = Model(Gurobi.Optimizer)
    if silent
        set_silent(model)
    end

    # apply T to all inputs
    df = df[T, :]

    ### GRID ###
    grid = Grid(network, connections, meta, T)
    add_grid_model(model=model,
        grid=grid,
        limit=limit,
        df=df,
        medians=medians,
        model_hp=model_hp
    )

    # sets 
    sets = get_sets(grid)
    B, H, H_HP, T = sets.B, sets.H, sets.H_HP, sets.T

    # variables
    P = model[:P]
    # P_HP = model[:P_HP]
    # PPD = model[:PPD]
    # J_heat = model[:J_c_heat]

    # define global model helper expressions
    @expressions(model, begin
        J_gen, sum(P[0, T]) # in Wh
        P_all_user_load[t in T], sum(P[i, t] for i in H)
    end)

    # define objective function
    # c_gen = 1E-2
    # @objective(model, Min, J_heat + c_gen * J_gen + J_ppd)
    @objective(model, Min, J_gen)

    # set solver options
    set_attribute(model, "BarHomogeneous", 1)
    set_attribute(model, "Presolve", 0)

    # optimize the model and return it
    optimize!(model)
    @assert is_solved_and_feasible(model)
    return model
end
