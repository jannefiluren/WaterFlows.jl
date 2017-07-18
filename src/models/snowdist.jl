
mutable struct SnowDistModel{s<:AbstractSnow} <: AbstractModel
    
    snow::s
    
end


function run_model(model::SnowDistModel, input::InputPT)

    nstep = size(input.prec, 2)

    swe_sim = zeros(nstep)

    for t in 1:nstep

        # Run snow model component

        set_input(model.snow, input, t)

        run_timestep(model.snow)

        swe_sim[t] = sum(model.snow.swe .* model.snow.frac)

    end

    return swe_sim

end


function set_input(snow::AbstractSnow, input::InputPT, t::Int64)

    for reg in eachindex(snow.frac)

        snow.tair[reg] = input.tair[reg, t]
        snow.p_in[reg] = input.prec[reg, t]

    end
    
end