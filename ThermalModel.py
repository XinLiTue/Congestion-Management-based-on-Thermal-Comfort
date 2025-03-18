def add_indoor_constraints(m, T_ind, Tem_ind, Heat, T_amb, C_house, R_house, PPD, lambdas, x_vals, y_vals, Pn, n_user, t):
    for i in range(n_user):
        # Indoor temperature dynamics
        if t == 0:
            m.addConstr(
                C_house * (T_ind[i, t] - Tem_ind[0]) == 1e3 * Heat[i, t] * 0.25 + (T_amb[0] - Tem_ind[0]) / R_house * 0.25,
                "IndoorTemChange0"
            )
        else:
            m.addConstr(
                C_house * (T_ind[i, t] - T_ind[i, t - 1]) == 1e3 * Heat[i, t] * 0.25 + (T_amb[t - 1] - T_ind[i, t - 1]) / R_house * 0.25,
                "IndoorTemChange"
            )

        # # Set indoor temperature
        # m.addConstr(T_ind[i, t] == Tem_ind[t], "IndoorTemSet")

        # PPD constraints
        # PPD(t)>= y_vals[j]+ (T_ind[i,t]-x_vals[j])*(y_vals[j+1]-y_vals[j])/(x_vals[j+1]-x_vals[j])
        for j in range(Pn - 1):
            m.addConstr(
                PPD[i, t] >= y_vals[j] + (T_ind[i, t] - x_vals[j]) * (y_vals[j + 1] - y_vals[j]) / (x_vals[j + 1] - x_vals[j]),
                f"PPD{j}"
            )
       
        