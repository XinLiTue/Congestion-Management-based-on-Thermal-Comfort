"""Constants used in the optimization problem."""

# p.u. base values
S_base = 1E5 # [VA]
V_base = 230.0 # [V]
Z_base = V_base^2 / S_base # [Ohm]
I_base = S_base / V_base # [A]

S_max = 1E6 * 0.23 * 0.16 * 3 / S_base # [pu]
const pf = 0.92 # baseload & HP power factor
const tan_phi_load = (sqrt(1 - pf^2)) / pf

# costs
const c_loss = 1.0
const c_grid = 500.0
const c_hp = 2000.0

# voltage constraints
const V_lb = 0.90 # [pu]
const V_ub = 1.10