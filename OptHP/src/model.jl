"""Define model related functions."""

function add_bus_variables(
    model::Model,
    grid::Grid
)
    # sets
    sets = get_sets(grid)
    B, H, H_HP, T = sets.B, sets.H, sets.H_HP, sets.T

    # variables
    @variables(model, begin
        # power injections
        P[T, B], (base_name = "PBusInjection")
        Q[T, B], (base_name = "QBusInjection")

        # voltage magnitude (3d)
        (V_lb * V_ref)^2 <= V[T, B] <= (V_ub * V_ref)^2, (base_name = "VoltSquare")

        # heat pump reactive power
        0 <= Q_HP[T, H_HP], (base_name = "HPreactivePower")
    end)
end


function add_line_variables(
    model::Model,
    grid::Grid
)
    # sets
    L = get_line_set(grid)  # all lines (start, end)
    T = get_sets(grid).T

    @variables(model, begin
        # power flows
        P_line[T, L], (base_name = "PLineFlow")
        Q_line[T, L], (base_name = "QLineFlow")

        # line current
        0 <= I_line[T, L], (base_name = "CurrentSquare")
    end)
end


function add_bus_constraints(
    model::Model,
    grid::Grid,
    limit::Vector{Float64}
)
    sets = get_sets(grid)
    B, H, H_HP, notH, T = sets.B, sets.H, sets.H_HP, sets.notH, sets.T

    # variables
    P = model[:P]
    Q = model[:Q]
    P_line = model[:P_line]
    Q_line = model[:Q_line]
    I_line = model[:I_line]

    # slack bus ID
    SB = get_slack_bus(grid).node

    # slack bus
    @constraints(model, begin
        [t in T], P[t, SB] >= limit[t], (base_name = "TransPowerLimitForCongestion")
        [t in T], P[t, SB]^2 + Q[t, SB]^2 <= S_trafo^2, (base_name = "TrafoLimit")
    end)

    for bus in grid.buses
        b = bus.node
        L_in = get_incoming_lines(b, grid)
        L_out = get_outgoing_lines(b, grid)

        # real power balance
        P_in = @expression(model, [t in T], sum(P_line[t, (l.start.node, l.stop.node)] for l in L_in))
        P_out = @expression(model, [t in T], sum(P_line[t, (l.start.node, l.stop.node)] for l in L_out))
        P_loss = @expression(model, [t in T], sum(I_line[t, (l.start.node, l.stop.node)] * l.R for l in L_in))

        # reactive power balance
        Q_in = @expression(model, [t in T], sum(Q_line[t, (l.start.node, l.stop.node)] for l in L_in))
        Q_out = @expression(model, [t in T], sum(Q_line[t, (l.start.node, l.stop.node)] for l in L_out))
        Q_loss = @expression(model, [t in T], sum(I_line[t, (l.start.node, l.stop.node)] * l.X for l in L_in))

        # if b=0 print everything
        if b == 0
            println("P_in: ", P_in)
            println("P_out: ", P_out)
            println("P_loss: ", P_loss)
            println("Q_in: ", Q_in)
            println("Q_out: ", Q_out)
            println("Q_loss: ", Q_loss)
            println("L_in: ", L_in)
            println("L_out: ", L_out)
        end

        @constraints(model, begin
            # real power balance (1a)
            [t in T], P[t, b] == P_in[t] - P_out[t] - P_loss[t], (base_name = "GridRealPowerBalance")
            # reactive power balance (1b)
            [t in T], Q[t, b] == Q_in[t] - Q_out[t] - Q_loss[t], (base_name = "GridReactivePowerBalance")
        end)
    end

    # load constraints for non-user nodes are set to zero
    @constraints(model, begin
        [t in T, h in notH; h ≠ 0], P[t, h] == 0, (base_name = "NonUserRealPowerBalance")
        [t in T, h in notH; h ≠ 0], Q[t, h] == 0, (base_name = "NonUserReactivePowerBalance")
    end)
end


function add_load_constraints(
    model::Model,
    grid::Grid,
    loads_real::DataFrame,
    loads_reactive::DataFrame
)
    sets = get_sets(grid)
    B, H, H_HP, notH, T = sets.B, sets.H, sets.H_HP, sets.notH, sets.T

    # variables
    P = model[:P]
    Q = model[:Q]
    P_base = loads_real
    Q_base = loads_reactive
    P_HP = model[:P_HP]
    Q_HP = model[:Q_HP]

    # constraints that apply *only* to non-heat pump users
    for bus in get_nonhp_buses(grid)
        b = bus.node
        @constraints(model, begin
            # load constraints (1a, 1b)
            [t in T], P[t, b] == P_base[t, "$b"], (base_name = "UserRealPowerBalance")
            [t in T], Q[t, b] == Q_base[t, "$b"], (base_name = "UserReactivePowerBalance")
        end)
    end

    # constraints that apply *only* to heat pump users
    for bus in get_hp_buses(grid)
        b = bus.node
        @constraints(model, begin
            # load constraints (1a, 1b)
            [t in T], P[t, b] == P_base[t, "$b"] + P_HP[t, b], (base_name = "P")
            [t in T], Q[t, b] == Q_base[t, "$b"] + Q_HP[t, b], (base_name = "Q")
            # heat pump constraints
            [t in T], Q_HP[t, b] == P_HP[t, b] * tan_phi_load, (base_name = "hp_pf")
        end)
    end
end


function add_line_constraints(
    model::Model,
    grid::Grid
)
    sets = get_sets(grid)
    B, H, H_HP, notH, T = sets.B, sets.H, sets.H_HP, sets.notH, sets.T
    L = get_line_set(grid)

    # variables
    P_line = model[:P_line]
    Q_line = model[:Q_line]
    I_line = model[:I_line]
    V = model[:V]

    for line in grid.lines
        l = (line.start.node, line.stop.node)
        @constraints(model, begin
            # voltage (1c)
            [t in T], V[t, l[1]] == V[t, l[2]] - 2 * (line.R * P_line[t, l] + line.X * Q_line[t, l]) +
                                    (line.R^2 + line.X^2) * I_line[t, l]
            # line current limit (1g)
            [t in T], I_line[t, l] <= line.Inom, (base_name = "LineCurrentLimit")
        end)
    end
end


function add_objectives(
    model::Model,
    grid::Grid
)
    sets = get_sets(grid)
    B, H, H_HP, notH, T = sets.B, sets.H, sets.H_HP, sets.notH, sets.T

    # variables
    I_line = model[:I_line]

    # heat pump energy cost 
    J_c = model[:J_c]

    # define objective functions
    @expressions(model, begin
        J_loss, sum(l.R * I_line[t, (l.start.node, l.stop.node)] for t in T, l in grid.lines)
        J_hp_c, sum(J_c[i] for i in H_HP)
    end)
end


function add_helper_expressions(model::Model, grid::Grid, loads_real::DataFrame)
    sets = get_sets(grid)
    B, H, H_HP, notH, T = sets.B, sets.H, sets.H_HP, sets.notH, sets.T

    # HP + baseload 
    @expression(model, P_load_user[t in T],
        sum(model[:P_HP][t, bus.node] for bus in get_hp_buses(grid)) +
        sum(loads_real[t, "$(bus.node)"] for bus in get_nonhp_buses(grid)))
end


function GEC(;
    network::DataFrame,
    connections::DataFrame,
    loads_real::DataFrame,
    loads_reactive::DataFrame,
    weather::DataFrame,
    limit::Tuple=(40:48, -15),
    meta::Dict=Dict(),
    silent=true,
    T=1:96
)
    # create the model
    model = Model(Gurobi.Optimizer)
    if silent
        set_silent(model)
    end

    # apply T to all dataframes 
    loads_real = loads_real[T, :]
    loads_reactive = loads_reactive[T, :]
    weather = weather[T, :]

    # build grid
    grid = build_grid(network, connections, meta, T)

    # setup congestion limit
    congestion_limit = ones(96) * (-50)
    congestion_limit[limit[1]] .= limit[2]

    ### HEAT PUMP ###
    add_heatpump_variables(model, grid, weather)
    add_heatpump_constraints(model, grid, weather)
    ### END HEAT PUMP ###

    ### GRID ###
    # grid variables
    add_bus_variables(model, grid)
    add_line_variables(model, grid)

    # grid constraints
    add_bus_constraints(model, grid, congestion_limit)
    add_line_constraints(model, grid)
    add_load_constraints(model, grid, loads_real, loads_reactive)
    ### END GRID ###


    # helper expressions
    add_helper_expressions(model, grid, loads_real)

    # objective functions
    add_objectives(model, grid)
    @objective(
        model,
        Min,
        model[:J_loss] * c_loss +
        model[:J_hp_c] * c_hp
    )

    # optimize the model and return it
    optimize!(model)
    @assert is_solved_and_feasible(model)
    return model
end
