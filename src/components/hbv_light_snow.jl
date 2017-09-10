

mutable struct HbvLightSnow
    
    state_SP::Array{Float64,2}
    state_WC::Array{Float64,2}

    par_TT::Array{Float64,1}
    par_CFMAX::Array{Float64,1}
    par_SFCF::Array{Float64,1}
    par_CFR::Array{Float64,1}
    par_CWH::Array{Float64,1}

    p_in::Array{Float64,1}
    tair::Array{Float64,1}
    q_out::Array{Float64,2}

    area::Array{Float64,2}
    
    tstep::Float64
    time::DateTime
    
end


function HbvLightSnow(tstep::Float64, time::DateTime, area::Array{Float64,2})
    
    @assert 1.0 <= tstep <= 24.0 "Time step outside allowed range (1.0 - 24.0h)"
    
    nveg, nelev = size(area)

    state_SP = zeros(nveg, nelev)
    state_WC = zeros(nveg, nelev)

    par_TT = fill(-1.0, nveg)
    par_CFMAX = fill(5.0, nveg)
    par_SFCF = fill(0.8, nveg)
    par_CFR = fill(0.05, nveg)
    par_CWH = fill(0.1, nveg)

    p_in = fill(0.0, nelev)
    tair = fill(0.0, nelev)
    q_out = fill(0.0, nveg, nelev)

    HbvLightSnow(state_SP, state_WC, par_TT, par_CFMAX, par_SFCF, 
    par_CFR, par_CWH, p_in, tair, q_out, area, tstep, time)

end


function get_param_ranges(m::HbvLightSnow)
    
    param_range = Dict(
    :par_TT => (-3.0, 3.0),
    :par_CFMAX => (1.0, 10.0),
    :par_SFCF => (0.5 ,2.0),
    :par_CFR => (0.0 ,0.1),
    :par_CWH => (0.0 ,0.2))
    
end


function init_states!(m::HbvLightSnow)
    
    m.state_SP .= zeros(m.state_SP)
    m.state_WC .= zeros(m.state_WC)
    
end


function run_timestep(m::HbvLightSnow)

    for elevzone in eachindex(m.p_in)

        P_zone = m.p_in[elevzone]
        T_zone = m.tair[elevzone]

        for vegzone = 1:size(m.area,1)

            SP = m.state_SP[vegzone, elevzone]
            WC = m.state_WC[vegzone, elevzone]

            tt    = m.par_TT[vegzone]
            CFMAX = m.par_CFMAX[vegzone]
            SFCF  = m.par_SFCF[vegzone]
            CFR   = m.par_CFR[vegzone]
            CWH   = m.par_CWH[vegzone]

            q_out = 0.0
            
            if SP > 0.0
                if P_zone > 0.0
                    if T_zone > tt
                        WC = WC + P_zone
                    else
                        SP = SP + P_zone * SFCF
                    end
                end # if P_zone
                if T_zone > tt
                    melt = CFMAX * (T_zone - tt)
                    if melt > SP
                        q_out = SP + WC
                        WC = 0.0
                        SP = 0.0
                    else
                        SP = SP - melt
                        WC = WC + melt
                        if WC >= CWH * SP
                            q_out = WC - CWH * SP
                            WC = CWH * SP
                        end
                    end
                else
                    refrez = CFR * CFMAX * (tt - T_zone)
                    if refrez > WC
                        refrez = WC
                    end
                    SP = SP + refrez
                    WC = WC - refrez
                end # if T_zone
            else
                if T_zone > tt
                    q_out = P_zone
                else
                    SP = P_zone * SFCF
                end # if T_zone
            end # if SP

            m.q_out[vegzone, elevzone] = q_out

        end

    end

    return nothing

end
