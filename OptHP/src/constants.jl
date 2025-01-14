"""Constants used in the optimization problem."""

const s_trafo = 0.23 * 0.16 * 3 # data from DACS area 7
const pf = 0.92 # baseload & HP power factor
const tan_phi_load = (sqrt(1 - pf^2)) / pf
const p_hp_max = 4.5 * 1e-3  # 4.5 [kW]
const p_hp_min = 1e-3 # 1 [kW]

# costs
const c_loss = 40
const c_grid = 500
const c_hp = 200
const c_pv = 400

# pv constraints
const pf_pv_limit = 0.95
const tan_phi_pv = (sqrt(1 - pf_pv_limit^2)) / pf_pv_limit

# voltage constraints
const V_ref = 0.23 # [kV]
const V_lb = 0.96 # [pu]
const V_ub = 1.04

# pv efficiency data for date "2024-02-01"
const hourly_data = [
    0, 0, 0, 0, 0, 0, 0, 0, 0.006, 0.053, 0.129, 0.179,
    0.166, 0.14, 0.094, 0.046, 0.007, 0, 0, 0, 0, 0, 0, 0
]
