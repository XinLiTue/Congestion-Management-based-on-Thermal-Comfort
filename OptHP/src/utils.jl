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
function convert_pu_base(model::Model, V_base::Float64, S_base::Float64)
    # get all variables
    vars = all_variables(model)

    # convert power units
    for (var, val) in vars
        if occursin("Power", var.base_name)
            val = val * S_base
        elseif occursin("Voltage", var.base_name)
            val = val * V_base
        end
    end
end