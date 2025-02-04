
# regression-based COP model 
COP(ΔT) = 20.3595 - 3.2061 * log2(1 + ΔT)


function add_heatpump_model(
    model::Model,
    grid::Grid,
    df::DataFrame;
    allow_boiler::Bool=true
)
    # sets
    sets = get_sets(grid)
    H_HP, T = sets.H_HP, sets.T
    HP = get_hp_buses(grid)
    Δt = grid.Δt

    # get electrical HP variable from grid model 
    # the variable has unit [p.u.], which we need to convert to kW
    P_HP = model[:P_HP] * S_base * 1E-3

    @variables(model, begin
        Φ_CV[H_HP, T], (base_name = "CVthermalPower")
        z_CV[H_HP, T], Bin, (base_name = "CVonoff")
        z_HP[H_HP, T], Bin, (base_name = "HPonoff")
        Φ_HP[H_HP, T], (base_name = "HPthermalPower")
        Te[H_HP, T, [:i, :e, :h]], (base_name = "StateSpaceModel")
    end)
    set_lower_bound.(Te[:, :, :h], 0.0)
    set_upper_bound.(Te[:, :, :h], 55)

    # if the boiler is not used, set z_CV to 0
    if !allow_boiler
        fix.(z_CV, 0; force=true)
    end

    # set initial values for the state variables
    for (k, v) in (
        :i => 20.0,
        :e => 20.0,
        :h => 20.0
    )
        fix.(Te[:, 1, k], v; force=true)
    end

    supply_T = Dict([bus.node => bus.T_supply for bus in HP])

    @expressions(model, begin
        # temperature delta between emission and ambient temperature
        ΔT[i in H_HP, t in T], supply_T[i] - df.T_a[t]
        η_COP[i in H_HP, t in T], COP(ΔT[i, t])

        # total thermal power output
        Φ_h[i in H_HP, t in T[1:end-1]], Φ_CV[i, t] + Φ_HP[i, t]

        # gas consumption s.t. Φ_CV[t] = g[t] * η_g * H_g / Δt
        g[i in H_HP, t in T], Φ_CV[i, t] * Δt / (η_g * H_g)

        # input vector
        u[i in H_HP, t in T[1:end-1]], [df.T_a[t], Φ_h[i, t], df.Φ_s[t]]

        # cost function for providing heat
        J_c_heat_t[i in H_HP, t in T], df.λ_e[t] * P_HP[i, t] * Δt + λ_g * g[i, t]
        J_c_heat[i in H_HP], sum(J_c_heat_t[i, t] for t in T)
    end)

    # add constraints for each heat pump
    for bus in HP
        i = bus.node
        A_d = bus.A
        B_d = bus.B

        @constraints(model, begin
            # operational constraints are in [kW]
            [t in T], Φ_CV[i, t] <= Φ_CV_max * z_CV[i, t]
            [t in T], Φ_CV[i, t] >= Φ_CV_min * z_CV[i, t]
            [t in T], Φ_HP[i, t] <= Φ_HP_max * z_HP[i, t]
            [t in T], Φ_HP[i, t] >= Φ_HP_min * z_HP[i, t]

            # COP 
            [t in T], Φ_HP[i, t] == η_COP[i, t] * P_HP[i, t]

            # comfort 
            [t in T], Te[i, t, :i] <= T_i_max
            [t in T], Te[i, t, :i] >= T_i_min

            # temperature at end can't be lower than at start 
            # Te[i, T[end], :i] >= Te[i, T[1], :i]

            # state-space model dynamics
            [t in T[1:end-1]], Te[i, t+1, :].data .== A_d * Te[i, t, :].data + B_d * u[i, t]
        end)
    end
end
