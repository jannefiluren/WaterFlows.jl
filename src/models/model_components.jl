

# Semidistributed model with snow and hydrological component

mutable struct ModelComp{s <: AbstractSnow, g <: AbstractGlacier, h <: AbstractSubsurf} <: AbstractModel
    
    snow::s
    glacier::g
    subsurf::h
    
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

@inline function set_input(s::AbstractSnow, input::InputPTE, t::Int64)

    for reg in 1:size(s.frac_lus, 2)
        s.tair[reg] = input.tair[reg, t]
        s.p_in[reg] = input.prec[reg, t]
    end
    
end


# Set input to glacier model

@inline function set_input(g::AbstractGlacier, input::InputPTE, t::Int64)

    for reg in eachindex(g.frac_lus)
        g.tair[reg] = input.tair[reg, t]
    end
    
end


@inline function set_input(g::NoGlacier, input::InputPTE, t::Int64)

    return nothing

end


# Set input to lumped subsurface model with glacier

@inline function set_input(h::AbstractSubsurfLumped, s::AbstractSnow, g::AbstractGlacier, input::InputPTE, t::Int64)

    h.epot = input.epot[t]
    
    h.p_in = 0.0

    for reg in eachindex(s.frac_lus)
        h.p_in += s.frac_lus[reg] * s.q_out[reg]
    end
    
    for reg in eachindex(g.frac_lus)
        h.p_in += g.frac_lus[reg] * g.q_out[reg]
    end
    
end


# Set input to lumped subsurface model without glacier

@inline function set_input(h::AbstractSubsurfLumped, s::AbstractSnow, g::NoGlacier, input::InputPTE, t::Int64)

    h.epot = input.epot[t]
    
    h.p_in = 0.0

    for reg in eachindex(s.frac_lus)
        h.p_in += s.frac_lus[reg] * s.q_out[reg]
    end
    
end


# Set input to distributed subsurface model with glacier

@inline function set_input(h::AbstractSubsurfDist, s::AbstractSnow, g::AbstractGlacier, input::InputPTE, t::Int64)

    h.epot = input.epot[t]

    nlus, nreg = size(h.frac_lus)

    h.p_in .= s.q_out

    for ireg = 1:nreg    # IS THIS CORRECT????
        h.p_in[g.iglacier, ireg] = g.q_out[ireg]
    end

    h.tth .= s.tth
    h.pcorr .= s.pcorr

    return nothing

end


# Set input to distributed subsurface model without glacier

@inline function set_input(h::AbstractSubsurfDist, s::AbstractSnow, g::NoGlacier, input::InputPTE, t::Int64)
    
    h.epot = input.epot[t]

    nlus, nreg = size(h.frac_lus)

    h.p_in  .= s.q_out
    h.tth   .= s.tth
    h.pcorr .= s.pcorr


    # input snow / no snow missing


    return nothing

end

