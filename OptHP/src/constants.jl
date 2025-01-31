"""Constants used in the optimization problem."""

### Electrical grid parameters ###

# p.u. base values
const S_base = 1E5 # [VA]
const V_base = 230.0 # [V]
const Z_base = V_base^2 / S_base # [Ohm]
const I_base = S_base / V_base # [A]

const S_max = 2 * 1E6 * 0.23 * 0.16 * 3 / S_base # [pu]
const pf = 0.92 # baseload & HP power factor
const tan_phi_load = (sqrt(1 - pf^2)) / pf

# costs
const c_loss = 1.0
const c_grid = 500.0
const c_hp = 2000.0

# voltage constraints
const V_lb = 0.90 # [pu]
const V_ub = 1.10

### Heat pump parameters ###

# parameters for costs
const λ_g = 1.286             # Gas price [€/m³] (https://www.energievergelijk.nl/energieprijzen/gasprijs)

# Operational constraints
const Φ_CV_min = 3.8            # Minimum thermal output power of the boiler [kW]
const Φ_CV_max = 14.8           # Maximum thermal output power of the boiler [kW]
const Φ_HP_min = 1.8            # Minimum thermal output power of the heat pump [kW]
const Φ_HP_max = 5.0            # Maximum thermal output power of the heat pump [kW]
const P_HP_max = 1.8            # Maximum electrical compressor power of the heat pump [kW]

# Efficiency parameters
const η_g = 0.95
const H_g = 9.77

# Comfort constraints
const T_i_max = 25.0           # Maximum upper comfort violation [°C]
const T_i_min = 18.0           # Maximum lower comfort violation [°C]