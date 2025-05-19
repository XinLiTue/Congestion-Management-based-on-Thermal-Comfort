def add_hhp_constraints(m, p_hp, h_hp, b_hp, g_boil, h_boil, b_boil, Heat, model_inf, n_user, t):
    
    for i in range(n_user):
        # Heat pump constraints
        m.addConstr(p_hp[i, t] <= b_hp[i, t] * model_inf.p_hp_max * model_inf.hp_own[i], f"hp_max{i,t}")
        m.addConstr(p_hp[i, t] >= b_hp[i, t] * model_inf.p_hp_min * model_inf.hp_own[i], f"hp_min{i,t}")
        m.addConstr(h_hp[i, t] == model_inf.COP * p_hp[i, t], f"hp_heat{i,t}")

        # Boiler constraints
        m.addConstr(h_boil[i, t] <= b_boil[i, t] * model_inf.p_boil_max, f"boil_max{i,t}")
        m.addConstr(h_boil[i, t] >= b_boil[i, t] * model_inf.p_boil_min, f"boil_min{i,t}")
        m.addConstr(1e3 * h_boil[i, t] == 4 * model_inf.gas_LHV * g_boil[i, t], f"boiler_heat{i,t}")

        # Switching constraint
        m.addConstr(b_boil[i, t] + b_hp[i, t] <= 1, f"HHP_constrain{i,t}")

        # Total heat output
        m.addConstr(Heat[i, t] == h_hp[i, t] + h_boil[i, t], f"HeatOutput{i,t}")