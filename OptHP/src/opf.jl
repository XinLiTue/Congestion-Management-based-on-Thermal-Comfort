
const s_trafo = 0.23 * 0.16 * 3 # data from DACS area 7
const pf = 0.92 # baseload & HP power factor
const tan_phi_load = (sqrt(1 - pf^2)) / pf
const p_hp_max = 4.5 * 1e-3  # 4.5 [kW]
const p_hp_min = 1e-3 # 1 [kW]

# costs
const c_loss = 40
const c_grid = 500
const c_hp_down = 200
const c_pv = 400

# pv constraints
const pf_pV_limit = 0.95
const tan_phi_pv = (sqrt(1 - pf_pV_limit^2)) / pf_pV_limit

# voltage constraints
const V_ref = 0.23 # [kV]
const V_lb = 0.96 # [pu]
const V_ub = 1.04

# pv efficiency data for date "2024-02-01"
const hourly_data = [
    0, 0, 0, 0, 0, 0, 0, 0, 0.006, 0.053, 0.129, 0.179,
    0.166, 0.14, 0.094, 0.046, 0.007, 0, 0, 0, 0, 0, 0, 0
]

function GEC(;
    network::DataFrame,
    connections::DataFrame,
    loads_real::DataFrame,
    loads_reactive::DataFrame,
    limit::Tuple=(40:48, -15 * 1e-3),
)
    # setup congestion limit
    congestion_limit = ones(96) * (-50 * 1e-3)
    congestion_limit[limit[1]] .= limit[2]

    # setup data
    pv_eff = repeat(hourly_data, inner=4)

    # check that loads / generation have equal length
    @assert nrow(loads_real) == nrow(loads_reactive) == length(pv_eff)

    # sets
    T = 1:length(pv_eff)                            # discrete time steps [-]
    B = 1:nrow(network)+1                           # all buses [-], slack bus is bus 1 (+1)
    H = connections[!, :Node]                       # user buses [-]
    notH = setdiff(B, H)                            # non-user buses [-]
    H_HP = connections[connections.HP.==1, :Node]   # buses with HP [-]

    # offset load data by user bus index
    P_base = OffsetArray(Array(loads_real), 1:length(T), maximum(notH)+1:maximum(H))
    Q_base = OffsetArray(Array(loads_reactive), 1:length(T), maximum(notH)+1:maximum(H))

    # PV capacity from DACS-HW data 
    pv_cap = connections[!, :PV] * 1e-3 # [kW]   
    pv_cap = OffsetArray(pv_cap, minimum(H):maximum(H))

    # create the model
    model = Model(Gurobi.Optimizer)
    set_silent(model)

    @variables(model, begin
        P[T, B], (base_name = "PBusInjection")
        Q[T, B], (base_name = "QBusInjection")

        # power flow end-1 because graph is a tree with |B| - 1 edges
        PLine[T, B[1:end-1]], (base_name = "PLine")
        QLine[T, B[1:end-1]], (base_name = "QLine")

        # squared voltage and current
        (V_lb * V_ref)^2 <= V[T, B] <= (V_ub * V_ref)^2, (base_name = "VoltSquare")
        0 <= I[T, B[1:end-1]], (base_name = "CurrentSquare")

        # photovoltaics
        0 <= P_pv[T, H], (base_name = "PVactivePower")
        Q_pv[T, H], (base_name = "PVreactivePower")
        0 <= P_pv_down[T, H], (base_name = "PVcurtailedActivePower")

        # heat pumps
        0 <= P_hp[T, H_HP], (base_name = "HPactivePower")
        0 <= Q_hp[T, H_HP], (base_name = "HPreactivePower")
        z_hp[T, H_HP], Bin, (base_name = "HPOnOff")
        0 <= P_hp_down[T, H_HP], (base_name = "HPcurtailedActivePower")

    end)

    # slack bus
    @constraints(model, begin
        [t in T], P[t, 1] >= congestion_limit[t], (base_name = "TransPowerLimitForCongestion")
        [t in T], P[t, 1]^2 + Q[t, 1]^2 <= s_trafo^2, (base_name = "TrafoLimit")
    end)

    # power balance of real power
    @constraint(model,
        [t in T, j in B],
        P[t, j] ==
        sum(PLine[t, i] for i in B[1:end-1] if network[i, :EndNode] == j) -                 # real power into bus
        sum(I[t, i] * network[i, :R] for i in B[1:end-1] if network[i, :EndNode] == j) -    # ohmic losses
        sum(PLine[t, i] for i in B[1:end-1] if network[i, :StartNode] == j)                 # power out of bus
    )

    # power balance of reactive power
    @constraint(model,
        [t in T, j in B],
        Q[t, j] ==
        sum(QLine[t, i] for i in B[1:end-1] if network[i, :EndNode] == j) -                 # reactive power into bus
        sum(I[t, i] * network[i, :X] for i in B[1:end-1] if network[i, :EndNode] == j) -    # reactive losses
        sum(QLine[t, i] for i in B[1:end-1] if network[i, :StartNode] == j)                 # reactive power out of bus
    )

    # electrical
    @constraints(model, begin

        # voltage relation
        [t in T, b in B[1:end-1]],
        V[t, network[b, :EndNode]] == V[t, network[b, :StartNode]] -
                                      2 * (network[b, :R] * PLine[t, b] + network[b, :X] * QLine[t, b]) +
                                      (network[b, :R]^2 + network[b, :X]^2) * I[t, b]

        # bus SOCP
        [t in T, b in B[1:end-1]],
        PLine[t, b]^2 + QLine[t, b]^2 <= V[t, network[b, :StartNode]] * I[t, b]

        # line current limit
        [t in T, b in B[1:end-1]], I[t, b] <= network[b, :Inom]

        # load constraints for non-user nodes are set to zero
        [t in T, h in notH[2:end]], P[t, h] == 0
        [t in T, h in notH[2:end]], Q[t, h] == 0

        # load constraints for users *without* a HP
        [t in T, h in H; h ∉ H_HP], P[t, h] == P_base[t, h] - P_pv[t, h]
        [t in T, h in H; h ∉ H_HP], Q[t, h] == Q_base[t, h] - Q_pv[t, h]

        # load constraints for users *with* a HP
        [t in T, h in H_HP], P[t, h] == P_base[t, h] - P_pv[t, h] + P_hp[t, h]
        [t in T, h in H_HP], Q[t, h] == Q_base[t, h] - Q_pv[t, h] + Q_hp[t, h]

        # photovoltaics
        [t in T, h in H], -P_pv[t, h] * tan_phi_pv <= Q_pv[t, h]    # pv_tan-
        [t in T, h in H], Q_pv[t, h] <= P_pv[t, h] * tan_phi_pv     # pv_tan+
        [t in T, h in H], P_pv[t, h] <= pv_cap[h] * pv_eff[t]       # pvMax
        [t in T, h in H], P_pv_down[t, h] == pv_cap[h] * pv_eff[t] - P_pv[t, h] # pv_down

        # heat pumps
        [t in T, h in H_HP], Q_hp[t, h] == P_hp[t, h] * tan_phi_load
        [t in T, h in H_HP], P_hp[t, h] <= z_hp[t, h] * p_hp_max # hp_max
        [t in T, h in H_HP], P_hp[t, h] >= z_hp[t, h] * p_hp_min # hp_min
        [t in T, h in H_HP], P_hp_down[t, h] >= p_hp_max - P_hp[t, h] # hp_down
    end)

    # define objective functions 
    @expressions(model, begin
        J_loss, sum(network[i, :R] * I[t, i]^2 for t in T, i in B[1:end-1])
        J_pv, sum(P_pv_down[t, i] for t in T, i in H)
        J_hp, sum(P_hp_down[t, i] for t in T, i in H_HP)
    end)

    @objective(model, Min, c_loss * J_loss + c_pv * J_pv + c_hp_down * J_hp)

    optimize!(model)
    @assert is_solved_and_feasible(model)

    return model
end