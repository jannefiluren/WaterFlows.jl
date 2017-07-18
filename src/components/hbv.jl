# Hydrological component

mutable struct Hbv <: AbstractHydro

    sm::Float64
    suz::Float64
    slz::Float64
    st_uh::Array{Float64,1}
    hbv_ord::Array{Float64,1}
    fc::Float64
    lp::Float64
    k0::Float64
    k1::Float64
    k2::Float64
    beta::Float64
    perc::Float64
    ulz::Float64
    maxbas::Float64
    p_in::Float64
    epot::Float64
    q_out::Float64
    tstep::Float64
    time::DateTime
    
end


function Hbv(tstep::Float64, time::DateTime)

    @assert 1.0 <= tstep <= 24.0 "Time step outside allowed range (1.0 - 24.0h)"

    sm     = 0.0
    suz    = 0.0
    slz    = 0.0
    st_uh  = zeros(Float64, 20)

    fc, lp, k0, k1, k2, beta, perc, ulz, maxbas = (100.0, 0.8, 0.05, 0.05, 0.01, 1.0, 2.0, 30., 2.5)

    p_in = 0.0
    epot = 0.0
    q_out = 0.0

    hbv_ord = compute_hbv_ord(maxbas)

    Hbv(sm, suz, slz, st_uh, hbv_ord, fc, lp, k0, k1, k2, beta, perc, ulz, maxbas, p_in, epot, q_out, tstep, time)

end


function get_param_ranges(model::Hbv)
    
    param_range = Dict(:fc => (1.0, 1000.0),
                    :lp => (0.5, 0.99),
                    :k0 => (0.001, 0.999),
                    :k1 => (0.001, 0.999),
                    :k2 => (0.001, 0.999),
                    :beta => (1.0, 5.0),
                    :perc => (0.1, 1000.0),
                    :ulz => (1.0, 1000.0),
                    :maxbas => (1.0, 20.0))
    
end










function init_states!(model::Gr4j)
    
    model.st[1] = 0.3 * model.x1
    model.st[2] = 0.5 * model.x3
    
    D = 1.30434782 * (model.tstep / 24.0) + 1.19565217
    
    UH1(model.ord_uh1, model.x4 * 24.0 / model.tstep, D)
    UH2(model.ord_uh2, model.x4 * 24.0 / model.tstep, D)  
    
    for i in eachindex(model.st_uh1)
        model.st_uh1[i] = 0.0
    end
    
    for i in eachindex(model.st_uh2)
        model.st_uh2[i] = 0.0
    end
    
end


function run_timestep(model::Gr4j)
    
    St     = model.st
    StUH1  = model.st_uh1
    StUH2  = model.st_uh2
    OrdUH1 = model.ord_uh1
    OrdUH2 = model.ord_uh2
    tstep  = model.tstep
    
    P1     = model.p_in
    E      = model.epot
    
    A = model.x1
    B = 0.9
    
    # Interception and production store
    
    if P1 <= E
        EN = E - P1
        PN = 0.0
        WS = EN / A
        if WS > 13.0
            WS = 13.0
        end
        TWS = tanh(WS)
        Sr = St[1] / A
        ER = St[1] * (2.0-Sr)*TWS/(1.0+(1.0-Sr)*TWS)
        AE = ER + P1
        St[1] = St[1] - ER
        PR = 0.0
    else
        EN = 0.0
        AE = E
        PN = P1 - E
        WS = PN / A
        if WS > 13.0
            WS = 13.0
        end
        TWS = tanh(WS)
        Sr = St[1] / A
        PS = A*(1.0-Sr*Sr)*TWS/(1.0+Sr*TWS)
        PR = PN - PS
        St[1] = St[1] + PS
    end
    
    # Percolation from production store
    
    if St[1] < 0.0
        St[1] = 0.0
    end
    
    scale_param = 25.62891 * (24.0 / tstep)
    
    Sr = St[1] / model.x1
    Sr = Sr * Sr
    Sr = Sr * Sr
    PERC = St[1] * (1.0-1.0/sqrt(sqrt(1.0 + Sr/scale_param)))
    
    St[1] = St[1] - PERC
    
    PR = PR + PERC
    
    # Split of effective rainfall into the two routing components
    
    PRHU1 = PR * B
    PRHU2 = PR * (1.0-B)
    
    # Convolution of unit hydrograph UH1
    
    NH = length(OrdUH1)
    
    for K = 1:NH-1
        StUH1[K] = StUH1[K+1] + OrdUH1[K]*PRHU1
    end
    
    StUH1[NH] = OrdUH1[NH] * PRHU1
    
    # Convolution of unit hydrograph UH2
    
    NH = length(OrdUH2)
    
    for K = 1:NH-1
        StUH2[K] = StUH2[K+1] + OrdUH2[K]*PRHU2
    end
    
    StUH2[NH] = OrdUH2[NH] * PRHU2
    
    # Potential intercatchment semi-exchange
    
    scale_param = (tstep / 24.0)
    
    Rr = St[2]/model.x3
    EXCH = scale_param * model.x2*Rr*Rr*Rr*sqrt(Rr)
    
    # Routing store
    
    AEXCH1 = EXCH
    
    if (St[2]+StUH1[1]+EXCH) < 0.0
        AEXCH1 = -St[2] - StUH1[1]
    end
    
    St[2] = St[2] + StUH1[1] + EXCH
    
    if St[2] < 0.0
        St[2] = 0.0
    end
    
    scale_param = (24.0 / tstep)
    Rr = St[2]^4 / (model.x3^4 * scale_param)
    
    QR = St[2] * (1.0-1.0/sqrt(sqrt(1.0+Rr)))
    
    St[2] = St[2] - QR
    
    # Runoff from direct branch QD
    
    AEXCH2 = EXCH
    
    if (StUH2[1]+EXCH) < 0.0
        AEXCH2 = -StUH2[1]
    end
    
    QD = max(0.0,StUH2[1]+EXCH)
    
    # Total runoff
    
    Q = QR + QD
    
    if Q < 0.0
        Q = 0.0
    end
    
    # Update struct
    
    model.time    += Dates.Hour(tstep)
    model.st      = St
    model.st_uh1  = StUH1
    model.st_uh2  = StUH2
    model.q_out   = Q
    
    return nothing
    
end





function compute_hbv_ord(maxbas)

  triang = Distributions.TriangularDist(0, maxbas)
  triang_cdf = Distributions.cdf(triang, 0:20)
  hbv_ord = diff(triang_cdf)

end

