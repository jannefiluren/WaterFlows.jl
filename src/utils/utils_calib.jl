

""" Single model run for calibration. """
function calib_wrapper(param, model::AbstractModel, input::AbstractInput, var_obs, warmup)

    @assert 1 <= warmup <= length(var_obs)

    set_params!(model, param)

    init_states!(model, input.time[1])

    var_sim = run_model(model, input)

    1.0 - nse(var_sim[warmup:end], var_obs[warmup:end])

end


""" Run model calibration. """
function run_model_calib(model::AbstractModel, input::AbstractInput, var_obs;
                         verbose = :silent, warmup = 3*365)

    param_range = get_param_ranges(model)

    calib_wrapper_tmp(param) = calib_wrapper(param, model, input, var_obs, warmup)

    res = bboptimize(calib_wrapper_tmp; SearchRange = param_range, TraceMode = verbose,  MaxSteps=5000)

    best_candidate(res)

end


""" Compute Nash-Sutcliffe efficiency. """
function nse(var_sim, var_obs)

    ikeep = .!isnan.(var_obs)
    
    var_sim = var_sim[ikeep]
    var_obs = var_obs[ikeep]
    
    1.0 .- sum((var_sim .- var_obs).^2) / sum((var_obs .- mean(var_obs)).^2)

end


""" Compute Kling-Gupta efficiency. """
function kge(var_sim, var_obs)
    
    if all(isnan, var_sim) || all(isnan, var_obs)

        kge = NaN

    else

        ikeep = .!isnan.(var_obs)

        var_sim = var_sim[ikeep]
        var_obs = var_obs[ikeep]

        r = cor(var_sim, var_obs)

        beta = mean(var_sim) / mean(var_obs)

        gamma = (std(var_sim) / mean(var_sim)) / (std(var_obs) / mean(var_obs))

        kge = 1 - sqrt( (r-1)^2 + (beta-1)^2 + (gamma-1)^2 )

    end

    return kge

end


""" Get parameter values for a model."""
function get_params(model::AbstractModel)
    res = []
    for name_comp in fieldnames(model)
        comp = getfield(model, name_comp)
        param_ranges = get_param_ranges(comp)
        for name_field in fieldnames(comp)
            if haskey(param_ranges, name_field)
                append!(res, getfield(comp, name_field))
            end
        end
    end
    return res
end


""" Get parameter ranges for a model."""
function get_param_ranges(model::AbstractModel)
    res = []
    for name_comp in fieldnames(typeof(model))
        comp = getfield(model, name_comp)
        param_ranges = get_param_ranges(comp)
        for name_field in fieldnames(typeof(comp))
            if haskey(param_ranges, name_field)
                param_tmp = getfield(comp, name_field)
                for i in 1:length(param_tmp)
                    push!(res, param_ranges[name_field])
                end
            end
        end
    end
    res = convert(Array{Tuple{Float64,Float64},1}, res)
    return res
end


"""Set parameters of a model."""
function set_params!(model::AbstractModel, param_input)
    
    # Check number of parameters    
    param_ranges = get_param_ranges(model)
    err_msg = "Model requires ($(length(param_ranges)) input parameters"
    @assert length(param_ranges) == length(param_input) 
    
    # Check parameter ranges
    for i in eachindex(param_input)
        param_min, param_max = param_ranges[i]
        err_msg = "Parameter number $i outside of allowed range $(param_ranges[i])"
        @assert param_min <= param_input[i] <= param_max err_msg
    end
    
    # Set parameter values
    iparam = 1
    for name_comp in fieldnames(typeof(model))
        comp = getfield(model, name_comp)
        param_ranges = get_param_ranges(comp)
        for name_field in fieldnames(typeof(comp))
            if haskey(param_ranges, name_field)
                param_tmp = getfield(comp, name_field)
                if param_tmp isa Float64
                    param_tmp = param_input[iparam]
                    iparam += 1
                else
                    for i in eachindex(param_tmp)
                        param_tmp[i] = param_input[iparam]
                        iparam += 1
                    end
                end
                setfield!(comp, name_field, param_tmp)
            end
        end
        setfield!(model, name_comp, comp)
    end
    return nothing
end


""" Initilize model state variables. """
function init_states!(model::AbstractModel, init_time::DateTime)

    for name_comp in fieldnames(typeof(model))

        comp = getfield(model, name_comp)

        init_states!(comp, init_time)

        setfield!(model, name_comp, comp)

    end
    
    return nothing

end

