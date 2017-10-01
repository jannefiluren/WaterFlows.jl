

# Semidistributed model with snow and hydrological component

mutable struct ModelComp{sn <: AbstractSnow, gl <: AbstractGlacier, ss <: AbstractSubsurf} <: AbstractModel
    
    snow::sn
    glacier::gl
    subsurf::ss
    
end

function run_model(model::ModelComp, input::InputPTE)

    nstep = size(input.prec, 2)

    q_out = zeros(nstep)

    for t in 1:nstep

        # Run snow model component

        set_input(model.snow, input, t)

        run_timestep(model.snow)

        # Run glacier model component

        set_input(model.glacier, input, t)

        run_timestep(model.glacier, model.snow)

        # Run subsurface model component

        set_input(model.subsurf, model.snow, model.glacier, input, t)

        run_timestep(model.subsurf)

        # Collect runoff as output

        q_out[t] = model.subsurf.q_out

    end

    return q_out

end


# Set input to snow model

@inline function set_input(sn::AbstractSnow, input::InputPTE, t::Int64)

    for reg in 1:size(sn.frac_lus, 2)
        sn.tair[reg] = input.tair[reg, t]
        sn.p_in[reg] = input.prec[reg, t]
    end
    
end


# Set input to glacier model

@inline function set_input(gl::AbstractGlacier, input::InputPTE, t::Int64)

    for reg in eachindex(gl.frac_lus)
        gl.tair[reg] = input.tair[reg, t]
    end
    
end


@inline function set_input(gl::NoGlacier, input::InputPTE, t::Int64)

    return nothing

end


# Set input to lumped subsurface model with glacier

@inline function set_input(ss::AbstractSubsurfLumped, sn::AbstractSnow, gl::AbstractGlacier, input::InputPTE, t::Int64)

    ss.epot = input.epot[t]
    
    ss.p_in = 0.0

    for reg in eachindex(sn.frac_lus)
        ss.p_in += sn.frac_lus[reg] * sn.q_out[reg]
    end
    
    for reg in eachindex(gl.frac_lus)
        ss.p_in += gl.frac_lus[reg] * gl.q_out[reg]
    end
    
end


# Set input to lumped subsurface model without glacier

@inline function set_input(ss::AbstractSubsurfLumped, sn::AbstractSnow, gl::NoGlacier, input::InputPTE, t::Int64)

    ss.epot = input.epot[t]
    
    ss.p_in = 0.0

    for reg in eachindex(sn.frac_lus)
        ss.p_in += sn.frac_lus[reg] * sn.q_out[reg]
    end
    
end


# Set input to distributed subsurface model with glacier

@inline function set_input(ss::AbstractSubsurfDist, sn::AbstractSnow, gl::AbstractGlacier, input::InputPTE, t::Int64)

    ss.epot = input.epot[t]
    
    for i in eachindex(ss.p_in)
        ss.p_in[i] = sn.q_out[i]
        ss.snow[i] = sn.swe[i] > 0.0
    end

    for i in eachindex(gl.q_out)
        ss.p_in[gl.iglacier, i] += gl.q_out[i]
    end

    return nothing

end


# Set input to distributed subsurface model without glacier

@inline function set_input(ss::AbstractSubsurfDist, sn::AbstractSnow, gl::NoGlacier, input::InputPTE, t::Int64)
    
    ss.epot = input.epot[t]

    for i in eachindex(ss.p_in)
        ss.p_in[i] = sn.q_out[i]
        ss.snow[i] = sn.swe[i] > 0.0
    end

    return nothing

end

