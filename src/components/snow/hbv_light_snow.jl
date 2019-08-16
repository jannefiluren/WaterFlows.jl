

mutable struct HbvLightSnow <: AbstractSnow
    
    swe::Array{Float64,2}
    whc::Array{Float64,2}
    tth::Array{Float64,1}
    ddf::Array{Float64,1}
    pcorr::Array{Float64,1}
    pfreeze::Array{Float64,1}
    pwhc::Array{Float64,1}
    p_in::Array{Float64,1}
    tair::Array{Float64,1}
    q_out::Array{Float64,2}
    frac_lus::Array{Float64,2}    
    tstep::Float64
    time::DateTime
    
end


function HbvLightSnow(tstep::Float64, time::DateTime, frac_lus::DataFrame)
    
    @assert tstep == 24.0 "Time step outside allowed range (24.0h)"
    
    frac_lus = convert(Array{Float64,2}, frac_lus)
    frac_lus = transpose(frac_lus)
    
    nlus, nreg = size(frac_lus)
    
    swe = zeros(nlus, nreg)
    whc = zeros(nlus, nreg)
    
    tth = fill(0.0, nlus)
    ddf = fill(5.0, nlus)
    pcorr = fill(1.0, nlus)
    pfreeze = fill(0.05, nlus)
    pwhc = fill(0.1, nlus)
    
    p_in = fill(0.0, nreg)
    tair = fill(0.0, nreg)
    q_out = fill(0.0, nlus, nreg)
    
    HbvLightSnow(swe, whc, tth, ddf, pcorr, 
    pfreeze, pwhc, p_in, tair, q_out, frac_lus, tstep, time)
    
end


function get_param_ranges(m::HbvLightSnow)
    
    param_range = Dict(
    :tth => (-3.0, 3.0),
    :ddf => (1.0, 10.0),
    :pcorr => (0.5 ,2.0),
    :pfreeze => (0.0 ,0.1),
    :pwhc => (0.0 ,0.2))
    
end


function init_states!(m::HbvLightSnow, init_time::DateTime)
    
    m.time = init_time
    
    m.swe .= zero(m.swe)
    m.whc .= zero(m.whc)
    
end


function get_water_stored(m::HbvLightSnow)

    water_stored = sum((m.swe + m.whc) .* m.frac_lus);

    return water_stored;

end


function run_timestep(m::HbvLightSnow)
    
    for ireg in eachindex(m.p_in)
        
        p_in = m.p_in[ireg]
        tair = m.tair[ireg]
        
        for ilus = 1:size(m.frac_lus, 1)
            
            if m.frac_lus[ilus, ireg] > 0.0
                
                swe = m.swe[ilus, ireg]
                whc = m.whc[ilus, ireg]
                
                tt      = m.tth[ilus]
                ddf     = m.ddf[ilus]
                pcorr   = m.pcorr[ilus]
                pfreeze = m.pfreeze[ilus]
                pwhc    = m.pwhc[ilus]
                
                q_out = 0.0
                
                if swe > 0.0
                    if p_in > 0.0
                        if tair > tt
                            whc = whc + p_in
                        else
                            swe = swe + p_in * pcorr
                        end
                    end
                    if tair > tt
                        melt = ddf * (tair - tt)
                        if melt > swe
                            q_out = swe + whc
                            whc = 0.0
                            swe = 0.0
                        else
                            swe = swe - melt
                            whc = whc + melt
                            if whc >= pwhc * swe
                                q_out = whc - pwhc * swe
                                whc = pwhc * swe
                            end
                        end
                    else
                        refrez = pfreeze * ddf * (tt - tair)
                        if refrez > whc
                            refrez = whc
                        end
                        swe = swe + refrez
                        whc = whc - refrez
                    end
                else
                    if tair > tt
                        q_out = p_in
                    else
                        swe = p_in * pcorr
                    end
                end
                
                m.swe[ilus, ireg] = swe
                m.whc[ilus, ireg] = whc
                
                m.q_out[ilus, ireg] = q_out
                
            end
            
        end
        
    end

    m.time += Dates.Hour(m.tstep)
    
    return nothing
    
end
