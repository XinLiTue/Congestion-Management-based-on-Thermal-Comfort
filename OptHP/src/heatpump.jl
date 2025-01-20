# Parameters for costs
const λ_g = 1.286             # Gas price [€/m³] (https://www.energievergelijk.nl/energieprijzen/gasprijs)
const c_u = 0.1               # Cost coefficient for upper comfort violation
const c_l = 10.0               # Cost coefficient for lower comfort violation

# Operational constraints
const Φ_CV_min = 3.8          # Minimum thermal output power of the boiler [kW]
const Φ_CV_max = 14.8         # Maximum thermal output power of the boiler [kW]
const Φ_HP_min = 1.8          # Minimum thermal output power of the heat pump [kW]
const Φ_HP_max = 5.0          # Maximum thermal output power of the heat pump [kW]
const P_HP_max = 1.8          # Maximum compressor power of the heat pump [kW]

# Efficiency parameters
const η_g = 0.95
const H_g = 9.77

# Comfort constraints
const T_i_max = 25.0           # Maximum upper comfort violation [°C]
const T_i_min = 20.0           # Maximum lower comfort violation [°C]


# default COP model 
COP(ΔT) = 20.3595 - 3.2061 * log2(1 + ΔT)


function add_heatpump_variables(
    model::Model,
    grid::Grid,
    df::DataFrame;
    use_boiler::Bool=true
)
    # sets
    sets = get_sets(grid)
    B, H, H_HP, notH, T = sets.B, sets.H, sets.H_HP, sets.notH, sets.T
    HP = get_hp_buses(grid)
    Δt = grid.Δt

    @variables(model, begin
        0 <= Φ_CV[T, H_HP] <= Φ_CV_max
        0 <= z_CV[T, H_HP] <= 1, Bin
        0 <= Φ_HP[T, H_HP] <= Φ_HP_max
        0 <= z_HP[T, H_HP] <= 1, Bin
        0 <= P_HP[T, H_HP] <= P_HP_max
        Te[T, H_HP, [:i, :e, :h]]   # indoor, envelope, emission temperature [°C]
    end)
    set_lower_bound.(Te[:, :, :h], 0.0)
    set_upper_bound.(Te[:, :, :h], 55)

    # if the boiler is not used, set z_CV to 0
    if !use_boiler
        fix(z_CV, 0; force=true)
    end

    # set initial values for the state variables
    for (k, v) in (
        :i => 21.0,
        :e => 21.0,
        :h => 21.0
    )
        fix.(Te[1, :, k], v; force=true)
    end

    supply_T = Dict([bus.node => bus.T_supply for bus in HP])

    @expressions(model, begin
        # temperature delta between emission and ambient temperature
        ΔT[t in T, b in H_HP], supply_T[b] - df.T_a[t]
        η_COP[t in T, b in H_HP], COP(ΔT[t, b])

        # total thermal power output
        Φ_h[t in T[1:end-1], b in H_HP], Φ_CV[t, b] + Φ_HP[t, b]

        # gas consumption s.t. Φ_CV[t] = g[t] * η_g * H_g / Δt
        g[t in T, b in H_HP], Φ_CV[t, b] * Δt / (η_g * H_g)

        # input vector
        u[t in T[1:end-1], b in H_HP], [df.T_a[t], Φ_h[t, b], df.Φ_s[t]]

        # cost functions
        J_c[b in H_HP], sum(df.λ_e[t] * P_HP[t, b] * Δt + λ_g * g[t, b] for t in T)
    end)
end


function add_heatpump_constraints(
    model::Model,
    grid::Grid,
    df::DataFrame
)
    # sets
    sets = get_sets(grid)
    B, H, H_HP, notH, T = sets.B, sets.H, sets.H_HP, sets.notH, sets.T
    HP = get_hp_buses(grid)

    # variables
    Φ_CV = model[:Φ_CV]
    Φ_HP = model[:Φ_HP]
    z_CV = model[:z_CV]
    z_HP = model[:z_HP]
    P_HP = model[:P_HP]
    η_COP = model[:η_COP]
    Te = model[:Te]
    u = model[:u]

    for bus in HP
        b = bus.node
        A_d = bus.A
        B_d = bus.B

        @constraints(model, begin
            # operational constraints
            [t in T], Φ_CV[t, b] <= Φ_CV_max * z_CV[t, b]
            [t in T], Φ_CV[t, b] >= Φ_CV_min * z_CV[t, b]
            [t in T], Φ_HP[t, b] <= Φ_HP_max * z_HP[t, b]
            [t in T], Φ_HP[t, b] >= Φ_HP_min * z_HP[t, b]

            # COP 
            [t in T], Φ_HP[t, b] == η_COP[t, b] * P_HP[t, b]

            # comfort 
            [t in T], Te[t, b, :i] <= T_i_max
            [t in T], Te[t, b, :i] >= T_i_min

            # state-space model dynamics
            [t in T[1:end-1]], Te[t+1, b, :].data .== A_d * Te[t, b, :].data + B_d * u[t, b]
        end)
    end
end