# Kavetski, D., and G. Kuczera (2007), Model smoothing strategies to remove microscale
# discontinuities and spurious secondary optima in objective functions in hydrological
# calibration, Water Resour. Res., 43, W03411, doi:10.1029/2006WR005195.
# 
# See Equation 8 for precipitation, and Equation 13 for melting/refreezing

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

@inline function pot_refreeze(tair, drf = 3.0, tth_refreeze = 0.0; m_refreeze = 0.5)
    
    t_m = -(tair - tth_refreeze) / m_refreeze
    pot_refreeze = drf * m_refreeze * (t_m + log(1. + exp(-t_m)))
    
    return pot_refreeze
    
end
