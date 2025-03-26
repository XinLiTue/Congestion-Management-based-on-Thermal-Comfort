def add_opf_constraints(m, PLine, QLine, l, v, p, q, model_inf, n_bus, t):
    # Slack bus constraints
    m.addQConstr(p[0, t] >= model_inf.congestion_limit[t], f"TransPowerLimitForCongestion{t}")
    m.addQConstr(p[0, t] * p[0, t] + q[0, t] * q[0, t] <= model_inf.s_trafo**2, f"trafoLimit{t}")

    # Non-slack bus power balance
    for j in range(n_bus):
        P_outFlow, P_inFlow, P_loss = 0, 0, 0
        for i in range(n_bus - 1):
            if model_inf.network.at[i, 'StartNode'] == j:
                P_outFlow += PLine[i, t]
            if model_inf.network.at[i, 'EndNode'] == j:
                P_inFlow += PLine[i, t]
                P_loss += l[i, t] * model_inf.network.at[i, 'R']
        m.addConstr(p[j, t] == -P_outFlow + (P_inFlow - P_loss), name=f"BusPower{j,t}")

    for j in range(n_bus):
        Q_outFlow, Q_inFlow, Q_loss = 0, 0, 0
        for i in range(n_bus - 1):
            if model_inf.network.at[i, 'StartNode'] == j:
                Q_outFlow += QLine[i, t]
            if model_inf.network.at[i, 'EndNode'] == j:
                Q_inFlow += QLine[i, t]
                Q_loss += l[i, t] * model_inf.network.at[i, 'X']
        m.addConstr(q[j, t] == -Q_outFlow + (Q_inFlow - Q_loss), name=f"BusReact{j,t}")

    # Voltage relation
    for i in range(n_bus - 1):
        m.addConstr(
            v[model_inf.network.at[i, 'EndNode'], t] == v[model_inf.network.at[i, 'StartNode'], t]
            - 2 * (model_inf.network.at[i, 'R'] * PLine[i, t] + model_inf.network.at[i, 'X'] * QLine[i, t])
            + (model_inf.network.at[i, 'R']**2 + model_inf.network.at[i, 'X']**2) * l[i, t],
            name=f"LineVoltage{i,t}"
        )

    # Second-order cone constraint
    for i in range(n_bus - 1):
        m.addQConstr(
            PLine[i, t] * PLine[i, t] + QLine[i, t] * QLine[i, t] <= l[i, t] * v[model_inf.network.at[i, 'StartNode'], t],
            name=f"BusSOCP{i,t}"
        )

    # Line current constraint
    for i in range(n_bus - 1):
        m.addConstr(l[i, t] <= (model_inf.network.at[i, 'Inom'])**2, name=f"LineCurrent{i,t}")

    # Slack bus voltage constraint
    m.addConstr(v[0, t] == model_inf.v_ref**2, f"transVol{t}")