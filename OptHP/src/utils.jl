const min_ms = Millisecond(Minute(1))


function interpolate_data(
    df::DataFrame,
    interval::Int64=1,
    t_col::Symbol=:timestamp;
    method::Symbol=:linear,
    include_interp::Bool=false
)
    # interval must be an integer larger than 1
    if interval < 1
        error("Interval must be an integer larger than 1")
    end

    # get minutes from timestamp
    minutes = Int64.(round.(cumsum(values(diff(df[!, t_col])) ./ min_ms), digits=0))
    minutes = vcat(0, minutes)

    if length(minutes) != nrow(df)
        println("Error: length of minutes is not equal to length of dataframe")
    end

    data_cols = Symbol.(filter(x -> x != String(t_col), names(df)))
    println("Interpolating columns: ", data_cols)
    data = Dict{Symbol,Vector}()

    # interpolate missing 
    if method == :linear
        itp_method = Linear()
    elseif method == :constant
        itp_method = Interpolations.Constant(Previous)
    else
        error("Interpolation method not recognized")
    end

    # define interpolation functions for each column
    for column in data_cols
        f_itp = interpolate((minutes,), df[!, column], Gridded(itp_method))
        data_inp = Float64.(f_itp(0:interval:minutes[end]))

        # constant columns are round to nearest integer (except T_room_set!)
        if method == :constant && column != :T_room_set
            data_inp = Int64.(round.(data_inp, digits=0))
        end
        data[column] = data_inp
    end

    # create new dataframe
    df_int = DataFrame(data)

    # interpolated timestamps
    df_int[!, t_col] = df[!, t_col][1] .+ Millisecond.(min_ms .* (0:interval:minutes[end]))

    # reorder columns
    df_int = df_int[:, [t_col; data_cols]]

    # indicate which rows are interpolated
    if include_interp
        interpolated = trues(nrow(df_int))
        interpolated[minutes.+1] .= false
        df_int.interpolated = interpolated
    end
    return df_int
end

# function to convert power units from p.u. to base units
function convert_pu_base(model::Model,
    V_base::Float64,
    S_base::Float64;
    slack_bus=55,
    user_bus=64,
    model_hp=false
)
    N = length(model[:P][1, :])
    J_heat = zeros(N)
    PPD = zeros(N)

    # convert variables to base units
    P_kw = Matrix(value.(model[:P]))[slack_bus, :] .* S_base .* 1E-3
    Q_kw = Matrix(value.(model[:Q]))[slack_bus, :] .* S_base .* 1E-3
    # J_heat
    if model_hp
        J_heat_euro = Vector(value.(model[:J_c_heat_t][user_bus, :]))
        PPD = Vector(value.(model[:PPD][user_bus, :]))
        T_i = Vector(value.(model[:Te])[user_bus, :, :i])
        T_e = Vector(value.(model[:Te])[user_bus, :, :e])
        T_h = Vector(value.(model[:Te])[user_bus, :, :h])
    end

    return DataFrame(
        P_kw=P_kw,
        Q_kw=Q_kw,
        J_heat_euro=J_heat_euro,
        PPD=PPD,
        T_i=T_i,
        T_e=T_e,
        T_h=T_h
    )
end

# save all results to a CSV file with the given filename
function save_result_csv(
    model::Model,
    path::String;
    model_hp::Bool=false,
)
    df = convert_pu_base(model, V_base, S_base, model_hp=model_hp)
    CSV.write(path, df)
end