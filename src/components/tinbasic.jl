# Snow component

mutable struct TinBasic <: AbstractSnow
    
    swe::Array{Float64,1}
    tth::Float64
    ddf::Float64
    pcorr::Float64
    p_in::Array{Float64,1}
    tair::Array{Float64,1}
    q_out::Array{Float64,1}
    frac::Array{Float64,1}
    tstep::Float64
    time::DateTime
    
end


function TinBasic(tstep::Float64, time::DateTime, frac::Array{Float64,1})
    
    swe   = zeros(Float64, length(frac))
    p_in  = zeros(Float64, length(frac))
    tair  = zeros(Float64, length(frac))
    q_out = zeros(Float64, length(frac))
    
    tth = 0.0
    ddf = 3.69
    pcorr = 1.02
    
    TinBasic(swe, tth, ddf, pcorr, p_in, tair, q_out, frac, tstep, time)
    
end


function get_param_ranges(model::TinBasic)
    
    param_range = Dict(:tth => (-3.0, 3.0),
                       :ddf => (0.1, 10.0),
                       :pcorr => (0.5, 2.0))
    
end


function init_states!(model::TinBasic)

  for ireg in eachindex(model.swe)
    model.swe[ireg] = 0.0
  end

end


function run_timestep(m::TinBasic)
    
    for reg in eachindex(m.swe)
        
        # Compute solid and liquid precipitation
        
        psolid, pliquid = split_prec(m.p_in[reg], m.tair[reg], m.tth)
        
        psolid  = psolid * m.pcorr
        
        # Compute snow melt
        
        M = pot_melt(m.tair[reg], m.ddf, m.tth)
        
        M = min(m.swe[reg],M)
        
        # Update snow water equivalents
        
        m.swe[reg] += psolid
        m.swe[reg] -= M
        
        # Compute snowpack runoff
        
        m.q_out[reg] = M + pliquid
        
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

