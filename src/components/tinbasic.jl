# Snow component

mutable struct TinBasic <: AbstractSnow
    
    swe::Array{Float64,2}
    tth::Float64
    ddf::Float64
    pcorr::Float64
    p_in::Array{Float64,1}
    tair::Array{Float64,1}
    q_out::Array{Float64,2}
    frac::Array{Float64,2}
    lus::Array{String,1}
    tstep::Float64
    time::DateTime
    
end


function TinBasic(tstep::Float64, time::DateTime, frac::Array{Float64,2}, lus::Array{String,1})
    
    swe   = zeros(Float64, size(frac))
    p_in  = zeros(Float64, size(frac,2))
    tair  = zeros(Float64, size(frac,2))
    q_out = zeros(Float64, size(frac))
    
    tth = 0.0
    ddf = 3.69
    pcorr = 1.02
    
    TinBasic(swe, tth, ddf, pcorr, p_in, tair, q_out, frac, lus, tstep, time)
    
end


function get_param_ranges(model::TinBasic)
    
    param_range = Dict(:tth => (-3.0, 3.0),
    :ddf => (0.1, 10.0),
    :pcorr => (0.5, 2.0))
    
end


function init_states!(model::TinBasic)
    
    for ilus in 1:size(model.frac, 1)
        for ireg in 1:size(model.frac, 2)
            model.swe[ilus, ireg] = 0.0
        end
    end
    
end


function run_timestep(m::TinBasic)
    
    for ireg in 1:size(m.frac, 2)

        # Compute solid and liquid precipitation
            
        psolid, pliquid = split_prec(m.p_in[ireg], m.tair[ireg], m.tth)
        
        psolid  = psolid * m.pcorr
        
        # Compute potential melt
            
        melt_pot = pot_melt(m.tair[ireg], m.ddf, m.tth)

        for ilus in 1:size(m.frac, 1)

            # Limit melt to available snow
            
            melt = min(m.swe[ilus, ireg], melt_pot)
            
            # Update snow water equivalents
            
            m.swe[ilus, ireg] += psolid - melt
            
            # Compute snowpack runoff
            
            m.q_out[ilus, ireg] = melt + pliquid
            
        end

    end
    
    m.time += Dates.Hour(m.tstep)
    
    return nothing
    
end

@inline function split_prec(prec, tair, tth_phase = 0.0; m_phase = 0.5)
    
    frac_snowfall = 1. / (1. + exp( (tair - tth_phase) / m_phase ))
    
    psolid = prec * frac_snowfall
    pliquid = prec - psolid
    
    return psolid, pliquid
    
end

@inline function pot_melt(tair, ddf = 3.0, tth_melt = 0.0; m_melt = 0.5)
    
    t_m = (tair - tth_melt) / m_melt
    pot_melt = ddf * m_melt * (t_m + log(1. + exp(-t_m)))
    
    return pot_melt
    
end

