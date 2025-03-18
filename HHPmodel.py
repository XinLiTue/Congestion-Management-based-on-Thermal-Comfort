def add_hhp_constraints(m, p_hp, q_hp, h_hp, b_hp, p_boil, h_boil, b_boil, Heat, gas_LHV, COP, n_user, t, p_hp_max, p_hp_min, p_boil_max, p_boil_min, hp_own):
    for i in range(n_user):
        # Heat pump constraints
        m.addConstr(p_hp[i, t] <= b_hp[i, t] * p_hp_max * hp_own[i], "hp_max")
        m.addConstr(p_hp[i, t] >= b_hp[i, t] * p_hp_min * hp_own[i], "hp_min")
        m.addConstr(h_hp[i, t] == COP * p_hp[i, t], "hp_heat")

        # Boiler constraints
        m.addConstr(h_boil[i, t] <= b_boil[i, t] * p_boil_max, "boil_max")
        m.addConstr(h_boil[i, t] >= b_boil[i, t] * p_boil_min, "boil_min")
        m.addConstr(1e3 * h_boil[i, t] == 4 * gas_LHV * p_boil[i, t], "boiler_heat")

        # Switching constraint
        m.addConstr(b_boil[i, t] + b_hp[i, t] <= 1, "HHP_constrain")

        # Total heat output
        m.addConstr(Heat[i, t] == h_hp[i, t] + h_boil[i, t], "HeatOutput")