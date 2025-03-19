import gurobipy as gp
from gurobipy import GRB,quicksum



def build_model(Time_day, model_inf, add_opf_constraints, add_hhp_constraints, add_indoor_constraints):
    

    

    
    n_bus=len(model_inf.network)+1
    user_index=model_inf.connect1["Node"].tolist()
    n_user= len(user_index)  #in DACS7, last one's load is 0

    pv_cap=model_inf.connect1["PV"]*1E-3 #pv capacity from dacs data
    pv_cap=pv_cap.tolist()

    hp_own=model_inf.connect1["HP"].tolist()

    m = gp.Model("GEC")
    m.Params.LogToConsole = 0

    # define variables
    p = m.addVars(n_bus,Time_day, lb = -float('inf'), vtype = GRB.CONTINUOUS, name = "PBusInjection")
    q = m.addVars(n_bus, Time_day,lb = -float('inf'),  vtype = GRB.CONTINUOUS, name = "QBusInjection")

    PLine = m.addVars(n_bus-1,Time_day, lb = -float('inf'),  vtype = GRB.CONTINUOUS, name = "PLine")
    QLine = m.addVars(n_bus-1,Time_day, lb = -float('inf'),  vtype = GRB.CONTINUOUS, name = "QLine")
    v = m.addVars(n_bus, Time_day,lb = (model_inf.v_lb*model_inf.v_ref)**2, ub = (model_inf.v_ub*model_inf.v_ref)**2, vtype = GRB.CONTINUOUS, name = "VoltSquare")
    l = m.addVars(n_bus-1, Time_day,lb = 0, vtype = GRB.CONTINUOUS, name = "CurrentSquare")
    p_pv = m.addVars(n_user,Time_day, lb = 0, vtype = GRB.CONTINUOUS, name = "PVactivePower")
    q_pv = m.addVars(n_user,Time_day, lb = -float('inf'),  vtype = GRB.CONTINUOUS, name = "PVreactivePower")

    p_hp = m.addVars(n_user,Time_day, lb = 0, vtype = GRB.CONTINUOUS, name = "HPactivePower")
    q_hp = m.addVars(n_user,Time_day, lb = 0,  vtype = GRB.CONTINUOUS, name = "HPreactivePower")

    p_pv_down = m.addVars(n_user,Time_day, lb = 0,  vtype = GRB.CONTINUOUS, name = "PVcurtailedActivePower")
    b_hp =m.addVars(n_user,Time_day,  vtype = GRB.BINARY, name = "HP_Open")
    p_hp_down = m.addVars(n_user,Time_day, lb = 0, vtype = GRB.CONTINUOUS, name = "HPcurtailedActivePower")

    #indoor variable
    h_hp=m.addVars(n_user,Time_day, lb = 0, vtype = GRB.CONTINUOUS, name = "HPHeat")
    g_boil = m.addVars(n_user,Time_day, lb = 0, vtype = GRB.CONTINUOUS, name = "GasToBoiler")
    h_boil = m.addVars(n_user,Time_day, lb = 0, vtype = GRB.CONTINUOUS, name = "BoilerHeat")
    b_boil =m.addVars(n_user,Time_day,  vtype = GRB.BINARY, name = "Boil_Open")

    Heat = m.addVars(n_user,Time_day, lb = 0, vtype = GRB.CONTINUOUS, name = "TotalHeat")
    T_ind = m.addVars(n_user,Time_day, lb = 0, vtype = GRB.CONTINUOUS, name = "IndoorTem")
    PPD  = m.addVars(n_user,Time_day, lb = 0, vtype = GRB.CONTINUOUS, name = "PPD")
    lambdas = m.addVars(model_inf.Pn, lb=0, ub=1, name="lambda")  # SOS2 variables

    # define constraints
    for t in range(Time_day):
            run_time=t
            pv_ef=model_inf.pvFactor.at[run_time,"Quarter_Hourly_Data"]
            p_baseload=model_inf.LoadPower.iloc[run_time,1:]*1E-3
            p_baseload=p_baseload.tolist()
            q_baseload=model_inf.LoadReact.iloc[run_time,1:]*1E-3
            q_baseload=q_baseload.tolist()

            add_opf_constraints(m, model_inf.network, PLine, QLine, l, v, p, q, n_bus, t, model_inf.congestion_limit, model_inf.s_trafo, model_inf.v_ref)
            add_hhp_constraints(m, p_hp, h_hp, b_hp, g_boil, h_boil, b_boil, Heat, model_inf.gas_LHV, model_inf.COP, n_user, t, model_inf.p_hp_max, model_inf.p_hp_min, model_inf.p_boil_max, model_inf.p_boil_min, hp_own)
            add_indoor_constraints(m, T_ind, model_inf.Tem_ind, Heat, model_inf.T_amb, model_inf.C_house, model_inf.R_house, PPD, lambdas, model_inf.x_vals, model_inf.y_vals, model_inf.Pn, n_user, t)

            for i in range(n_bus):
                if i not in user_index and i != 0:
                    m.addConstr(p[i,t] == 0,"busP=0")
                    m.addConstr(q[i,t] == 0,"busQ=0")
                if i in user_index:
                    m.addConstr(p[i,t] ==  p_baseload[i-(n_bus-n_user)] - p_pv[i-(n_bus-n_user),t] + p_hp[i-(n_bus-n_user),t], "LoadP")
                    m.addConstr(q[i,t] ==  q_baseload[i-(n_bus-n_user)]  - q_pv[i-(n_bus-n_user),t] + q_hp[i-(n_bus-n_user),t], "LoadQ")

            for i in range(n_user):        
                m.addConstr(q_pv[i,t] == p_pv[i,t] * model_inf.tan_phi_pv,"pv_tan")
                m.addConstr(q_hp[i,t] == p_hp[i,t] * model_inf.tan_phi_load,"hp_tan")         
                m.addConstr(p_pv[i,t] == pv_cap[i]*pv_ef,"pvMax")

    # define objective
    obj = 0\
    +quicksum(1e3*model_inf.ele_price[int(t/4)]*(-p[0,t])*0.25  for t in range(Time_day) )\
    + quicksum(model_inf.gas_price[int(t/4)]*g_boil[i,t] for t in range(Time_day) for i in range(n_user))\
    + quicksum(model_inf.PPD_Price*PPD[i,t] for t in range(Time_day) for i in range(n_user))\
    # + quicksum(l[i,t]*network.at[i,'R'] * 1e3*ele_price[int(t/4)]*0.25 for t in range(Time_day) for i in range(n_bus-1))\
    # + quicksum(c_hp_down*p_hp_down[i,t]*0.25 for t in range(Time_day) for i in range(n_user))\
    # + quicksum(c_pv*p_pv_down[i,t]*0.25 for t in range(Time_day) for i in range(n_user))\


  
    m.setObjective(obj, GRB.MINIMIZE)

    m.Params.MIPGap = 0.005
    m.optimize()


    print(f"Optimization Runtime: {m.Runtime} seconds")
    if m.status == GRB.OPTIMAL:
        print("Model solved successfully!")    
        # m.write("opf_model.sol")
    else:
        print(f"Model status: {m.status}")


    dict_optimizedResults = {
    "p": [[p[i,t].X  for t in range(Time_day)]for i in range(n_bus)],
    "q": [[q[i,t].X  for t in range(Time_day)] for i in range(n_bus)], 
    "v_value": [[v[i,t].X  for t in range(Time_day)] for i in range(n_bus)],   
    "p_hp_down": [[p_hp_down[i,t].X for t in range(Time_day)]for i in range(n_user)], 
    "p_hp": [[p_hp[i,t].X for t in range(Time_day)]for i in range(n_user)],
    "p_pv_down": [[p_pv_down[i,t].X for t in range(Time_day)]for i in range(n_user)], 
    "p_pv": [[p_pv[i,t].X for t in range(Time_day)]for i in range(n_user)],
    "T_ind": [[T_ind[i,t].X for t in range(Time_day)]for i in range(n_user)],
    "g_boil": [[g_boil[i,t].X for t in range(Time_day)]for i in range(n_user)],
    "PPD": [[PPD[i,t].X for t in range(Time_day)]for i in range(n_user)]
    } 


    return m,dict_optimizedResults