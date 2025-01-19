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
        t_u[T, H_HP] >= 0
        t_l[T, H_HP] >= 0
        Te[T, H_HP, [:i, :e, :h]]   # indoor, envelope, emission temperature [°C]
        η_COP[T, H_HP] >= 1.0
        ΔT[T, H_HP]
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
        # fix(T[1, k], v; force=true)
        fix.(Te[1, :, k], v; force=true)
    end

    @expressions(model, begin
        # temperature delta between emission and ambient temperature
        Φ_h[t in T[1:end-1], b in H_HP], Φ_CV[t, b] + Φ_HP[t, b]

        # gas consumption s.t. Φ_CV[t] = g[t] * η_g * H_g / Δt
        g[t in T, b in H_HP], Φ_CV[t, b] * Δt / (η_g * H_g)

        # input vector
        u[t in T[1:end-1], b in H_HP], [df.T_a[t], Φ_h[t, b], df.Φ_s[t]]

        # # cost functions
        # J_c, sum(df.λ_e[t] * P_HP[t] * Δt + λ_g * g[t] for t in T)
        # J_d, c_u * sum(t_u) + c_l * sum(t_l)
    end)
end