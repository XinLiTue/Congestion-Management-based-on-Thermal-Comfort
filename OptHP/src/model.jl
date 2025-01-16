"""Define model related functions."""

function add_bus_variables(
    model::Model,
    grid::Grid,
    T::UnitRange{Int},
)
    B = Set(getfield.(grid.buses, :node))
    H = Set([bus.node for bus in grid.buses if bus isa User])
    H_HP = Set([bus.node for bus in grid.buses if bus isa HeatPumpBus])

    return @variables(model, begin
        # power injections
        P[T, B], (base_name = "PBusInjection")
        Q[T, B], (base_name = "QBusInjection")

        # voltage magnitudes
        (V_lb * V_ref)^2 <= V[T, B] <= (V_ub * V_ref)^2, (base_name = "VoltSquare")

        # photovoltaics
        0 <= P_pv[T, H], (base_name = "PVactivePower")
        Q_pv[T, H], (base_name = "PVreactivePower")
        0 <= P_pv_down[T, H], (base_name = "PVcurtailedActivePower")

        # heat pumps (TEST, simple model)
        0 <= P_hp[T, H_HP], (base_name = "HPactivePower")
        0 <= Q_hp[T, H_HP], (base_name = "HPreactivePower")
        z_hp[T, H_HP], Bin, (base_name = "HPOnOff")
        0 <= P_hp_down[T, H_HP], (base_name = "HPcurtailedActivePower")
    end)
end


function GEC(;
    network::DataFrame,
    connections::DataFrame,
    loads_real::DataFrame,
    loads_reactive::DataFrame,
    limit::Tuple=(40:48, -15 * 1e-3),
    meta::Dict=Dict()
)


    # create the model
    model = Model(Gurobi.Optimizer)
    set_silent(model)

    # setup data
    pv_eff = repeat(hourly_data, inner=4)

    # discrete time steps [-]
    T = 1:length(pv_eff)

    # build grid
    grid = build_grid(network, connections, meta)

    # add variables to the model
    add_bus_variables(model, grid, T)

end
