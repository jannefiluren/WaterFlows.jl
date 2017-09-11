
mutable struct HbvLightSubsurf <: AbstractSubsurfDist
    
    state_SM::Array{Float64,2}
    SUZ::Float64
    SLZ::Float64
    state_uh::Array{Float64,1}
    hbv_ord::Array{Float64,1}
    
    PERC::Float64
    K0::Float64
    K1::Float64
    K2::Float64
    UZL::Float64
    par_TT::Array{Float64,1}
    par_SFCF::Array{Float64,1}
    par_FC::Array{Float64,1}
    par_LP::Array{Float64,1}
    par_BETA::Array{Float64,1}
    par_MAXBAS::Float64
    
    snow::Array{Bool,2}

    p_in::Array{Float64,2}
    tair::Array{Float64,2}
    epot::Float64
    q_out::Float64
    
    lake::Float64                 # This needs to be handled!!!
    frac::Array{Float64,2}
    
    tstep::Float64
    time::DateTime
    
end



function HbvLightSubsurf(tstep::Float64, time::DateTime, frac_lus::DataFrame, lake::Float64)
    
    @assert 1.0 <= tstep <= 24.0 "Time step outside allowed range (1.0 - 24.0h)"

    @assert lake == 0.0 "No lakes allowed currently"
    
    frac = convert(Array{Float64,2}, frac_lus)
    frac = transpose(frac)

    nveg, nelev = size(frac)
    
    state_SM = zeros(nveg, nelev)
    SUZ = 0.0
    SLZ = 0.0
    
    PERC = 0.7
    K0 = 0.2
    K1 = 0.08
    K2 = 0.03
    UZL = 20.0
    
    par_FC = fill(250.0, nveg)
    par_LP = fill(0.7, nveg)
    par_BETA = fill(3.0, nveg)
    par_TT = fill(0.0, nveg)      # Temperature threshold (set equal to snow model in connector function)
    par_SFCF = fill(1.0, nveg)     # Snowfall correction factor (set equal to snow model in connector function)
    par_MAXBAS = 2.5
    
    snow = fill(false, nveg, nelev)

    p_in = fill(0.0, nveg, nelev)
    tair = fill(0.0, nveg, nelev)
    epot = 0.0
    q_out = 0.0
    
    hbv_ord = compute_hbv_ord(par_MAXBAS)
    state_uh = zeros(hbv_ord)
    
    HbvLightSubsurf(state_SM, SUZ, SLZ, state_uh, hbv_ord, PERC, K0, K1, K2, UZL, par_TT, par_SFCF,
    par_FC, par_LP, par_BETA, par_MAXBAS, snow, p_in, tair, epot, q_out, lake, frac, tstep, time)
    
end



function get_param_ranges(m::HbvLightSubsurf)
    
    param_range = Dict(
    :PERC => (0.1, 1000.0),
    :K0 => (0.001, 0.999),
    :K1 => (0.001, 0.999),
    :K2 => (0.001, 0.999),
    :ULZ => (1.0, 1000.0),
    :par_FC => (50.0, 500.0),
    :par_LP => (0.3, 1.0),
    :par_BETA => (1.0, 6.0),
    :par_MAXBAS => (1.0, 7.0))
    
end

function init_states!(m::HbvLightSubsurf)
    
    m.state_SM .= zeros(m.state_SM)
    m.SUZ = 0.0
    m.SLZ = 0.0
    
    m.hbv_ord = compute_hbv_ord(m.par_MAXBAS)
    
    m.state_uh .= 0.0
    
end


function run_timestep(m::HbvLightSubsurf)
    
    pET = m.epot
    
    pot_E = pET * (1 - m.lake)
    
    till_Qsum = 0.0
    avg_act_E = 0.0
    
    for elevzone = 1:size(m.frac, 2)
        
        for vegzone = 1:size(m.frac, 1)
            
            if m.frac[vegzone, elevzone] > 0.0
                
                SM = m.state_SM[vegzone, elevzone]
               
                FC   = m.par_FC[vegzone]
                LP   = m.par_LP[vegzone]
                BETA = m.par_BETA[vegzone]
                
                insoil = m.p_in[vegzone, elevzone]
                
                #soil routine
                
                till_Q = 0.0
                old_SM = SM
                if insoil > 0.0
                    if insoil < 1.0
                        y = insoil
                    else
                        mi = floor(insoil)   # IS THIS CORRECT?
                        y = insoil - mi
                        for i in 1:mi
                            dQdP = (SM / FC) ^ BETA
                            if dQdP > 1.0
                                dQdP = 1.0
                            end
                            SM = SM + 1.0 - dQdP
                            till_Q = till_Q + dQdP
                        end
                    end
                    dQdP = (SM / FC) ^ BETA
                    if dQdP > 1.0
                        dQdP = 1.0
                    end
                    SM = SM + (1 - dQdP) * y
                    till_Q = till_Q + dQdP * y
                end # if insoil
                
                mean_SM = (SM + old_SM) / 2.0
                if mean_SM < (LP * FC)
                    act_E = pot_E * mean_SM / (LP * FC)
                else
                    act_E = pot_E
                end
                #if snow                       # Currently snow does not influence actual evapotranspiration
                #    act_E = 0.0
                #end
                SM = SM - act_E
                if SM < 0.0
                    SM = 0.0
                end
                
                avg_act_E = avg_act_E + act_E * m.frac[vegzone, elevzone]
                till_Qsum = till_Qsum + till_Q * m.frac[vegzone, elevzone]
                m.state_SM[vegzone, elevzone] = SM
                
            end # if frac
            
        end # veg loop
        
    end # elev loop
    
    till_box1 = till_Qsum
    
    # generation of runoff
    m.SUZ = m.SUZ + till_box1
    if ( m.SUZ - m.PERC / (1 - m.lake) ) < 0.0
        m.SLZ = m.SLZ + m.SUZ * (1.0 - m.lake)
        m.SUZ = 0.0
    else
        m.SLZ = m.SLZ + m.PERC
        m.SUZ = m.SUZ - m.PERC / (1.0 - m.lake)
    end
    
    tt = mean(m.par_TT)   # CHECK WITH JAN WHAT THIS MEANS
    SFCF = mean(m.par_SFCF)   # CHECK WITH JAN WHAT THIS MEANS
    
    if mean(m.tair) > tt   # CHECK WITH JAN WHAT THIS MEANS
        m.SLZ = max(m.SLZ - pET * m.lake, 0.0)
        avg_act_E = avg_act_E + min(m.SLZ, pET * m.lake)
    end
    if mean(m.tair) <= tt   # CHECK WITH JAN WHAT THIS MEANS
        m.SLZ = m.SLZ + SFCF * mean(m.p_in) * m.lake   # CHECK WITH JAN WHAT THIS MEANS
    else
        m.SLZ = m.SLZ + mean(m.p_in) * m.lake   # CHECK WITH JAN WHAT THIS MEANS
    end
    
    Q_box1 = m.K1 * m.SUZ
    if m.SUZ < m.UZL
        Q_box0 = 0.0
    else
        Q_box0 = m.K0 * (m.SUZ - m.UZL)
    end
    
    Q_box2 = m.K2 * m.SLZ
    m.SUZ = m.SUZ - Q_box1 - Q_box0
    m.SLZ = m.SLZ - Q_box2
    Q_gen = Q_box1 + Q_box2 + Q_box0
    
    # Convolution of unit hydrograph
    
    nh = length(m.hbv_ord)
    for k = 1:nh-1
        m.state_uh[k] = m.state_uh[k+1] + m.hbv_ord[k] * Q_gen
    end
    m.state_uh[nh] = m.hbv_ord[nh] * Q_gen
    
    # Output runoff
    
    m.q_out = m.state_uh[1]
    
    return nothing    
    
end    


