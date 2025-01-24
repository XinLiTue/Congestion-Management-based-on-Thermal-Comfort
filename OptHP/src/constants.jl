"""Constants used in the optimization problem."""

const S_max = 1E3 * 0.23 * 0.16 * 3 # [kVA]
const pf = 0.92 # baseload & HP power factor
const tan_phi_load = (sqrt(1 - pf^2)) / pf

# costs
const c_loss = 1.0
const c_grid = 500.0
const c_hp = 2000.0
const c_pv = 400.0

# voltage constraints
const V_ref = 0.23 # [kV]
const V_lb = 0.96 # [pu]
const V_ub = 1.04

# pv efficiency data for date "2024-02-01"
const hourly_data = [
    0, 0, 0, 0, 0, 0, 0, 0, 0.006, 0.053, 0.129, 0.179,
    0.166, 0.14, 0.094, 0.046, 0.007, 0, 0, 0, 0, 0, 0, 0
]
