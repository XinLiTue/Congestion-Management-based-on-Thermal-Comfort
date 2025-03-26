from gurobipy import GRB, quicksum
import numpy as np
### State-Space Model Dynamics

# The indoor temperature evolution is modeled using a state-space formulation:

# $$
# \mathbf{T}_{i,t+1} = \mathbf{A} \mathbf{T}_{i,t} + \mathbf{B} \mathbf{U}_{i,t}
# $$

# where:

# - $ \mathbf{T}_{i,t} $ is the state vector at time $ t $:

# $$
# \mathbf{T}_{i,t} = 
# \begin{bmatrix} 
# T_{i,t} \\ 
# T_{h,t} \\ 
# T_{a,t} 
# \end{bmatrix}
# $$

# - $ \mathbf{U}_{i,t} $ is the input vector:

# $$
# \mathbf{U}_{i,t} = 
# \begin{bmatrix} 
# T_a \\ 
# Heat_{i,t} \\ 
# Solar[t]
# \end{bmatrix}
# $$

# - $ \mathbf{A} $ is the state transition matrix:

# $$
# \mathbf{A} =
# \begin{bmatrix} 
# A_{11} & A_{12} & A_{13} \\
# A_{21} & A_{22} & A_{23} \\
# A_{31} & A_{32} & A_{33} 
# \end{bmatrix}
# $$

# - $ \mathbf{B} $ is the input matrix:

# $$
# \mathbf{B} =
# \begin{bmatrix} 
# B_{11} & B_{12} & B_{13} \\
# B_{21} & B_{22} & B_{23} \\
# B_{31} & B_{32} & B_{33} 
# \end{bmatrix}
# $$

# - $ \mathbf{T}_{i,t} = \begin{bmatrix} T_i \\ T_h \\ T_a \end{bmatrix}_t $ represents the state variables:
#   - $ T_i $: Interior temperature
#   - $ T_h $: Emission system temperature
#   - $ T_a $: House envelope temperature, assume to equal to ambient temperature



def add_indoor_constraints(m, T_ind, model_inf, Heat, PPD, n_user, t, IFRC=True):
    for i in range(n_user):
        # Indoor temperature dynamics    

        A=[[0.9754234372338041,0.013301970222414965,0.004636748139716821],
            [0.0034206080570929574,0.9964772560231249,8.382306645063422e-6],
            [0.20432236921074795,0.0014364106308368987,0.793524864333977]]
 

        B=[[0.006637844404064138,0.0022459832989580805,0.16643854382482537],    
        [9.37536131372134e-5,2.6505141147745418e-6,0.006442730239255753],
        [0.0007163558244381842,0.8289496956732726,0.017961087687400356]]

        # T_i is T_ind[i,t], T_a is T_amb[t]
        T_h=25
        Solar = np.ones(96)
        if IFRC:
        #simple RC model
            if t == 0:
                m.addConstr(
                    model_inf.C_house * (T_ind[i, t] - model_inf.Tem_ind[0]) == 
                    1e3 * Heat[i, t] * 0.25 + (model_inf.T_amb[0] - model_inf.Tem_ind[0]) / model_inf.R_house * 0.25,
                    f"IndoorTemChange0{i}"
                )
            else:
                m.addConstr(
                    model_inf.C_house * (T_ind[i, t] - T_ind[i, t - 1]) == 
                    1e3 * Heat[i, t] * 0.25 + (model_inf.T_amb[t - 1] - T_ind[i, t - 1]) / model_inf.R_house * 0.25,
                    f"IndoorTemChange{i,t}"
                )
        else:
            # State-space model using A B U T
            if t == 0:
                m.addConstr(
                    T_ind[i, t] == A[0][0] * model_inf.Tem_ind[0] + A[0][1] * T_h + A[0][2] * model_inf.T_amb[t] + 
                                  B[0][0] * model_inf.T_amb[t] + B[0][1] * 1e3* Heat[i, t] + B[0][2] * model_inf.solar_output[t],
                    f"IndoorTemChange0{i}"
                )
            else:
                m.addConstr(
                    T_ind[i, t] == A[0][0] * T_ind[i, t - 1] + A[0][1] * T_h + A[0][2] * model_inf.T_amb[t] + 
                                  B[0][0] * model_inf.T_amb[t] + B[0][1] * 1e3* Heat[i, t] + B[0][2] * model_inf.solar_output[t],
                    f"IndoorTemChange{i,t}"
                )




        # # Set indoor temperature
        # m.addConstr(T_ind[i, t] == Tem_ind[t], "IndoorTemSet")

        # PPD constraints
        # PPD(t)>= y_vals[j]+ (T_ind[i,t]-x_vals[j])*(y_vals[j+1]-y_vals[j])/(x_vals[j+1]-x_vals[j])
        # SOS2 constraints
        # m.addConstr(sum(lambdas[j] for j in range(Pn)) == 1, "SumToOne")
        # m.addConstr(T_ind[i, t] == sum(x_vals[j] * lambdas[j] for j in range(Pn)), "T_indEq")
        # m.addConstr(PPD[i, t] == sum(y_vals[j] * lambdas[j] for j in range(Pn)), "PPDEq")
        # m.addSOS(GRB.SOS_TYPE2, [lambdas[j] for j in range(Pn)])
        Pn = len(model_inf.x_vals)
        for j in range(Pn - 1):
            m.addConstr(
                PPD[i, t] >= model_inf.y_vals[j] + (T_ind[i, t] - model_inf.x_vals[j]) * 
                                      (model_inf.y_vals[j + 1] - model_inf.y_vals[j]) / 
                                      (model_inf.x_vals[j + 1] - model_inf.x_vals[j]),
                f"PPD_Conic{j,t}"
            )
       
        