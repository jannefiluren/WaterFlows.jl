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