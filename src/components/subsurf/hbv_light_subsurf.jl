
mutable struct HbvLightSubsurf <: AbstractSubsurfDist
    
    sm::Array{Float64,2}
    suz::Float64
    slz::Float64
    st_uh::Array{Float64,1}
    ord_uh::Array{Float64,1}    
    perc::Float64
    k0::Float64
    k1::Float64
    k2::Float64
    uzl::Float64
    fc::Array{Float64,1}
    lp::Array{Float64,1}
    beta::Array{Float64,1}
    maxbas::Float64
    snow::Array{Bool,2}
    p_in::Array{Float64,2}
    epot::Float64
    q_out::Float64
    frac_lus::Array{Float64,2}
    tstep::Float64
    time::DateTime
    
end



function HbvLightSubsurf(tstep::Float64, time::DateTime, frac_lus::DataFrame)
    
    frac_lus = convert(Array{Float64,2}, frac_lus)
    frac_lus = transpose(frac_lus)

    nlus, nreg = size(frac_lus)
    
    sm = zeros(nlus, nreg)
    suz = 0.0
    slz = 0.0
    
    perc = 0.7
    k0 = 0.2
    k1 = 0.08
    k2 = 0.03
    uzl = 20.0
    
    fc = fill(250.0, nlus)
    lp = fill(0.7, nlus)
    beta = fill(3.0, nlus)
    maxbas = 2.5
    
    snow = fill(false, nlus, nreg)

    p_in = fill(0.0, nlus, nreg)
    epot = 0.0
    q_out = 0.0
    
    ord_uh = compute_hbv_ord(maxbas)
    st_uh = zero(ord_uh)
    
    HbvLightSubsurf(sm, suz, slz, st_uh, ord_uh, perc, k0, k1, k2, uzl, fc,
    lp, beta, maxbas, snow, p_in, epot, q_out, frac_lus, tstep, time)
    
end



function get_param_ranges(m::HbvLightSubsurf)
    
    param_range = Dict(
    :perc => (0.1, 1000.0),
    :k0 => (0.001, 0.999),
    :k1 => (0.001, 0.999),
    :k2 => (0.001, 0.999),
    :ulz => (1.0, 1000.0),
    :fc => (50.0, 500.0),
    :lp => (0.3, 1.0),
    :beta => (1.0, 6.0),
    :maxbas => (1.0, 7.0))
    
end


function init_states!(m::HbvLightSubsurf, init_time::DateTime)
    
    m.time = init_time
    
    m.sm .= zero(m.sm)
    m.suz = 0.0
    m.slz = 0.0
    
    m.ord_uh = compute_hbv_ord(m.maxbas)
    
    m.st_uh .= 0.0
    
end


function run_timestep(m::HbvLightSubsurf)
    
    epot = m.epot
    
    to_qsum = 0.0
    avg_aet = 0.0
    
    for ireg = 1:size(m.frac_lus, 2)
        
        for ilus = 1:size(m.frac_lus, 1)
            
            if m.frac_lus[ilus, ireg] > 0.0
                
                sm   = m.sm[ilus, ireg]
                fc   = m.fc[ilus]
                lp   = m.lp[ilus]
                beta = m.beta[ilus]
                
                insoil = m.p_in[ilus, ireg]
                
                to_q = 0.0
                old_sm = sm
                
                if insoil > 0.0
                    if insoil < 1.0
                        y = insoil
                    else
                        mi = floor(insoil)   # IS THIS CORRECT?
                        y = insoil - mi
                        for i in 1:mi
                            dqdp = (sm / fc) ^ beta
                            if dqdp > 1.0
                                dqdp = 1.0
                            end
                            sm = sm + 1.0 - dqdp
                            to_q = to_q + dqdp
                        end
                    end
                    dqdp = (sm / fc) ^ beta
                    if dqdp > 1.0
                        dqdp = 1.0
                    end
                    sm = sm + (1 - dqdp) * y
                    to_q = to_q + dqdp * y
                end
                
                mean_sm = (sm + old_sm) / 2.0
                if mean_sm < (lp * fc)
                    aet = epot * mean_sm / (lp * fc)
                else
                    aet = epot
                end
                #if snow                       # Currently snow does not influence actual evapotranspiration
                #    aet = 0.0
                #end
                sm = sm - aet
                if sm < 0.0
                    sm = 0.0
                end
                
                avg_aet = avg_aet + aet * m.frac_lus[ilus, ireg]
                to_qsum = to_qsum + to_q * m.frac_lus[ilus, ireg]
                m.sm[ilus, ireg] = sm
                
            end
            
        end
        
    end
    
    to_suz = to_qsum
    
    # generation of runoff
    m.suz = m.suz + to_suz
    if ( m.suz - m.perc ) < 0.0
        m.slz = m.slz + m.suz
        m.suz = 0.0
    else
        m.slz = m.slz + m.perc
        m.suz = m.suz - m.perc
    end
           
    q_box1 = m.k1 * m.suz
    if m.suz < m.uzl
        q_box0 = 0.0
    else
        q_box0 = m.k0 * (m.suz - m.uzl)
    end
    
    q_box2 = m.k2 * m.slz
    m.suz = m.suz - q_box1 - q_box0
    m.slz = m.slz - q_box2
    q_gen = q_box1 + q_box2 + q_box0
    
    # Convolution of unit hydrograph
    
    nh = length(m.ord_uh)
    for k = 1:nh-1
        m.st_uh[k] = m.st_uh[k+1] + m.ord_uh[k] * q_gen
    end
    m.st_uh[nh] = m.ord_uh[nh] * q_gen
    
    # Output runoff
    
    m.q_out = m.st_uh[1]

    # Update time

    m.time += Dates.Hour(m.tstep)
    
    return nothing
    
end    


