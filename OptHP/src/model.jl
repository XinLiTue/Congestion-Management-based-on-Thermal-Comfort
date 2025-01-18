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

        # photovoltaics
        0 <= P_pv[T, H], (base_name = "PVactivePower")
        Q_pv[T, H], (base_name = "PVreactivePower")
        0 <= P_pv_down[T, H], (base_name = "PVcurtailedActivePower")

        # heat pumps (TEST, simple model)
        0 <= P_hp[T, H_HP], (base_name = "HPactivePower")
        0 <= Q_hp[T, H_HP], (base_name = "HPreactivePower")
        z_hp[T, H_HP], Bin, (base_name = "HPOnOff")
        0 <= P_hp_down[T, H_HP], (base_name = "HPcurtailedActivePower")
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

        @constraints(model, begin
            # real power balance (3a)
            [t in T], P[t, b] == P_in[t] - P_out[t] - P_loss[t], (base_name = "GridRealPowerBalance")
            # reactive power balance (3b)
            [t in T], Q[t, b] == Q_in[t] - Q_out[t] - Q_loss[t], (base_name = "GridReactivePowerBalance")
        end)
    end

    # load constraints for non-user nodes are set to zero
    @constraints(model, begin
        [t in T, h in notH; h ≠ 0], P[t, h] == 0
        [t in T, h in notH; h ≠ 0], Q[t, h] == 0
    end)
end


function add_load_constraints(
    model::Model,
    grid::Grid,
    loads_real::DataFrame,
    loads_reactive::DataFrame,
    pv_eff::Vector{Float64}
)
    sets = get_sets(grid)
    B, H, H_HP, notH, T = sets.B, sets.H, sets.H_HP, sets.notH, sets.T

    # variables
    P = model[:P]
    Q = model[:Q]
    P_base = loads_real
    Q_base = loads_reactive
    P_pv = model[:P_pv]
    Q_pv = model[:Q_pv]
    P_hp = model[:P_hp]
    Q_hp = model[:Q_hp]
    z_hp = model[:z_hp]
    P_pv_down = model[:P_pv_down]
    P_hp_down = model[:P_hp_down]

    # constraints that apply *only* to non-heat pump users
    for bus in get_nonhp_buses(grid)
        b = bus.node
        @constraints(model, begin
            # load constraints
            [t in T], P[t, b] == P_base[t, "$b"] - P_pv[t, b], (base_name = "UserRealPowerBalance")
            [t in T], Q[t, b] == Q_base[t, "$b"] - Q_pv[t, b], (base_name = "UserReactivePowerBalance")
        end)
    end

    # constraints that apply *only* to heat pump users
    for bus in get_hp_buses(grid)
        b = bus.node
        @constraints(model, begin
            # load constraints 
            [t in T], P[t, b] == P_base[t, "$b"] - P_pv[t, b] + P_hp[t, b]
            [t in T], Q[t, b] == Q_base[t, "$b"] - Q_pv[t, b] + Q_hp[t, b]
            # heat pump constraints
            [t in T], Q_hp[t, b] == P_hp[t, b] * tan_phi_load, (base_name = "hp_pf")
            [t in T], P_hp[t, b] <= z_hp[t, b] * P_hp_max, (base_name = "hp_max")
            [t in T], P_hp[t, b] >= z_hp[t, b] * p_hp_min, (base_name = "hp_min")
            [t in T], P_hp_down[t, b] >= P_hp_max - P_hp[t, b], (base_name = "hp_down")
        end)

        ### TEST CONSTRAINT ###
        # SET P_hp_down = 0
        @constraints(model, begin
            [t in T], P_hp_down[t, b] == 10.0  # basically a dummy constraint
        end)
    end

    # constraints that apply to *all* users
    for bus in get_user_buses(grid)
        b = bus.node
        PV_cap = bus.PV
        @constraints(model, begin
            [t in T], -P_pv[t, b] * PV_tan_phi <= Q_pv[t, b], (base_name = "pv_pf-")
            [t in T], Q_pv[t, b] <= P_pv[t, b] * PV_tan_phi, (base_name = "pv_pf+")
            [t in T], P_pv[t, b] <= PV_cap * pv_eff[t], (base_name = "pv_max")
            [t in T], P_pv_down[t, b] == PV_cap * pv_eff[t] - P_pv[t, b], (base_name = "pv_down")
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
            # voltage (3c)
            [t in T], V[t, l[1]] == V[t, l[2]] - 2 * (line.R * P_line[t, l] + line.X * Q_line[t, l]) +
                                    (line.R^2 + line.X^2) * I_line[t, l]
            # line current limit (3g)
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
    P_pv_down = model[:P_pv_down]
    P_hp_down = model[:P_hp_down]
    I_line = model[:I_line]

    # define objective functions
    @expressions(model, begin
        J_loss, sum(l.R * I_line[t, (l.start.node, l.stop.node)] for t in T, l in grid.lines)
        J_pv, sum(P_pv_down[t, i] for t in T, i in H)
        J_hp, sum(P_hp_down[t, i] for t in T, i in H_HP)
    end)
end


function add_helper_expressions(model::Model, grid::Grid, loads_real::DataFrame)
    sets = get_sets(grid)
    B, H, H_HP, notH, T = sets.B, sets.H, sets.H_HP, sets.notH, sets.T

    # HP + baseload 
    @expression(model, P_load_user[t in T],
        sum(model[:P_hp][t, bus.node] for bus in get_hp_buses(grid)) +
        sum(loads_real[t, "$(bus.node)"] for bus in get_nonhp_buses(grid)))
end


function GEC(;
    network::DataFrame,
    connections::DataFrame,
    loads_real::DataFrame,
    loads_reactive::DataFrame,
    limit::Tuple=(40:48, -15 * 1e-3),
    meta::Dict=Dict(),
    silent=true
)
    # create the model
    model = Model(Gurobi.Optimizer)
    if silent
        set_silent(model)
    end

    # setup data
    pv_eff = repeat(hourly_data, inner=4)

    # discrete time steps [-]
    T = 1:length(pv_eff)

    # build grid
    grid = build_grid(network, connections, meta, T)

    # setup congestion limit
    congestion_limit = ones(96) * (-50 * 1e-3)
    congestion_limit[limit[1]] .= limit[2]

    # variables
    add_bus_variables(model, grid)
    add_line_variables(model, grid)

    # constraints
    add_bus_constraints(model, grid, congestion_limit)
    add_line_constraints(model, grid)
    add_load_constraints(model, grid, loads_real, loads_reactive, pv_eff)

    # objective functions
    add_objectives(model, grid)
    @objective(model, Min, model[:J_loss] * c_loss +
                           model[:J_pv] * c_pv +
                           model[:J_hp] * c_hp
    )

    # helper expressions
    add_helper_expressions(model, grid, loads_real)

    # optimize the model and return it
    optimize!(model)
    @assert is_solved_and_feasible(model)
    return model
end
