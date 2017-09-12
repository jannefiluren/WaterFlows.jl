

# Semidistributed model with snow and hydrological component

mutable struct SemiDistComp{s <: AbstractSnow, g <: AbstractGlacier, h <: AbstractSubsurf} <: AbstractModel
    
    snow::s
    glacier::g
    hydro::h
    
end

function run_model(model::SemiDistComp, input::InputPTE)

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

        set_input(model.hydro, model.snow, model.glacier, input, t)

        run_timestep(model.hydro)

        # Collect runoff as output

        q_out[t] = model.hydro.q_out

    end

    return q_out

end


# Set input to snow model

function set_input(s::AbstractSnow, input::InputPTE, t::Int64)

    for reg in 1:size(s.frac, 2)
        s.tair[reg] = input.tair[reg, t]
        s.p_in[reg] = input.prec[reg, t]
    end
    
end


# Set input to glacier model

function set_input(g::AbstractGlacier, input::InputPTE, t::Int64)

    for reg in eachindex(g.frac)
        g.tair[reg] = input.tair[reg, t]
    end
    
end


# Set input to lumped subsurface model with glacier

function set_input(h::AbstractSubsurfLumped, s::AbstractSnow, g::AbstractGlacier, input::InputPTE, t::Int64)

    h.epot = input.epot[t]
    
    h.p_in = 0.0

    for reg in eachindex(s.frac)
        h.p_in += s.frac[reg] * s.q_out[reg]
    end
    
    for reg in eachindex(g.frac)
        h.p_in += g.frac[reg] * g.q_out[reg]
    end
    
end


# Set input to lumped subsurface model without glacier

function set_input(h::AbstractSubsurfLumped, s::AbstractSnow, g::NoGlacier, input::InputPTE, t::Int64)

    h.epot = input.epot[t]
    
    h.p_in = 0.0

    for reg in eachindex(s.frac)
        h.p_in += s.frac[reg] * s.q_out[reg]
    end
    
end


# Set input to distributed subsurface model with glacier

function set_input(h::AbstractSubsurfDist, s::AbstractSnow, g::AbstractGlacier, input::InputPTE, t::Int64)

    h.epot = input.epot[t]

    nlus, nreg = size(h.frac)

    h.p_in .= 0.0

    for ilus = 1:nlus, ireg = 1:nreg
        h.p_in[ilus, ireg] += s.frac[ilus, ireg] * s.q_out[ilus, ireg]
    end

    for ireg = 1:nreg
        h.p_in[g.iglacier, ireg] = g.frac[ireg] * g.q_out[ireg]
    end

    h.tth .= s.tth
    h.pcorr .= s.pcorr

    return nothing

end


# Set input to distributed subsurface model without glacier

function set_input(h::AbstractSubsurfDist, s::AbstractSnow, g::NoGlacier, input::InputPTE, t::Int64)
    
    h.epot = input.epot[t]

    nlus, nreg = size(h.frac)

    h.p_in .= 0.0

    for ilus = 1:nlus, ireg = 1:nreg
        h.p_in[ilus, ireg] += s.frac[ilus, ireg] * s.q_out[ilus, ireg]
    end

    h.tth .= s.tth
    h.pcorr .= s.pcorr


    # input snow / no snow missing


    return nothing

end








#= 
# Semidistributed model with full model

mutable struct SemiDistFull{h<:HbvLight} <: AbstractModel
    
    hydro::h
    
end

function run_model(model::SemiDistFull, input::InputPTE)

    nstep = size(input.prec, 2)

    q_out = zeros(nstep)

    for t in 1:nstep

        # Run full model

        set_input(model.hydro, input, t)

        run_timestep(model.hydro)

        q_out[t] = model.hydro.q_out

    end

    return q_out

end

function set_input(model::HbvLight, input::InputPTE, t::Int64)

    for reg in 1:size(model.area,2)

        model.tair[reg] = input.tair[reg, t]
        model.p_in[reg] = input.prec[reg, t]

    end
    
    model.epot = input.epot[t]

end =#
