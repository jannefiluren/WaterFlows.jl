

""" Single model run for calibration. """
function calib_wrapper(param, model::AbstractModel, input::AbstractInput, q_obs, warmup)

    @assert 1 <= warmup <= length(q_obs)

    set_params!(model, param)

    init_states!(model)

    q_sim = run_model(model, input)

    1.0 - nse(q_sim[warmup:end], q_obs[warmup:end])

end


""" Run model calibration. """
function run_model_calib(model::AbstractModel, input::AbstractInput, q_obs;
                         verbose = :silent, warmup = 3*365)

    param_range = get_param_ranges(model)

    calib_wrapper_tmp(param) = calib_wrapper(param, model, input, q_obs, warmup)

    res = bboptimize(calib_wrapper_tmp; SearchRange = param_range, TraceMode = verbose)

    best_candidate(res)

end


""" Compute Nash-Sutcliffe efficiency. """
function nse(q_sim, q_obs)

    ikeep = .!isnan.(q_obs)
    
    q_sim = q_sim[ikeep]
    q_obs = q_obs[ikeep]
    
    1.0 - sum((q_sim - q_obs).^2) / sum((q_obs - mean(q_obs)).^2)

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
    for name_comp in fieldnames(model)
        comp = getfield(model, name_comp)
        param_ranges = get_param_ranges(comp)
        for name_field in fieldnames(comp)
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
    for name_comp in fieldnames(model)
        comp = getfield(model, name_comp)
        param_ranges = get_param_ranges(comp)
        for name_field in fieldnames(comp)
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
function init_states!(model::AbstractModel)

    for name_comp in fieldnames(model)

        comp = getfield(model, name_comp)

        init_states!(comp)

        setfield!(model, name_comp, comp)

    end
    
    return nothing

end

