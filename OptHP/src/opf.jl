
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
    limit::Tuple=(40:48, -15 * 1e-3),
)
    # setup congestion limit
    congestion_limit = ones(96) * (-50 * 1e-3)
    congestion_limit[limit[1]] .= limit[2]

    # setup data
    expanded_data = repeat(hourly_data, inner=4)

    # sets
    H = connections[!, :Node]       # user nodes [1, 2, 3, ...]
    T = 1:length(expanded_data)     # discrete time steps [1, 2, 3, ...]
    B = 1:nrow(network)+1           # bus nodes [1, 2, 3, ...]

    # create the model
    model = Model(Gurobi.Optimizer)
    set_silent(model)

    @variables(model, begin
        P[B, T], (base_name = "PBusInjection")
        Q[B, T], (base_name = "QBusInjection")

        # power flow
        PLine[B[1:end-1], T], (base_name = "PLine")
        QLine[B[1:end-1], T], (base_name = "QLine")

        # squared voltage and current
        (V_lb * V_ref)^2 <= V[B, T] <= (V_ub * V_ref)^2, (base_name = "VoltSquare")
        0 <= I[B[1:end-1], T], (base_name = "CurrentSquare")

        # photovoltaics
        0 <= P_pv[H, T], (base_name = "PVactivePower")
        Q_pv[H, T], (base_name = "PVreactivePower")

        # heat pumps
        0 <= P_hp[H, T], (base_name = "HPactivePower")
        0 <= Q_hp[H, T], (base_name = "HPreactivePower")

        # curtailed power
        0 <= P_pv_down[H, T], (base_name = "PVcurtailedActivePower")
        b_hp[H, T], Bin, (base_name = "HP_Open")
        0 <= P_hp_down[H, T], (base_name = "HPcurtailedActivePower")
    end)

    # slack bus
    @constraints(model, begin
        [t in T], P[1, t] >= congestion_limit[t], (base_name = "TransPowerLimitForCongestion")
        [t in T], P[1, t]^2 + Q[1, t]^2 <= s_trafo^2, (base_name = "TrafoLimit")
    end)

    # non-slack bus, power balance
end