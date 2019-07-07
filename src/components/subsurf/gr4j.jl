# Hydrological component

mutable struct Gr4j <: AbstractSubsurfLumped
    
    st::Array{Float64,1}
    st_uh1::Array{Float64,1}
    st_uh2::Array{Float64,1}
    ord_uh1::Array{Float64,1}
    ord_uh2::Array{Float64,1}
    x1::Float64
    x2::Float64
    x3::Float64
    x4::Float64
    p_in::Float64
    epot::Float64
    q_out::Float64
    exch::Float64
    tstep::Float64
    time::DateTime
    
end


function Gr4j(tstep::Float64, time::DateTime)
    
    @assert 1.0 <= tstep <= 24.0 "Time step outside allowed range (1.0 - 24.0h)"
    
    n_ord = ceil(Int64, 20.0 * 24.0 / tstep)
    
    st      = zeros(Float64, 2)
    st_uh1  = zeros(Float64, n_ord)
    st_uh2  = zeros(Float64, 2 * n_ord)
    ord_uh1 = zeros(Float64, n_ord)
    ord_uh2 = zeros(Float64, 2 * n_ord)
    
    x1, x2, x3, x4  = (74.59, 0.81, 214.98, 1.24)
    
    st[1] = 0.3 * x1
    st[2] = 0.5 * x3
    
    D = 1.30434782 * (tstep / 24.0) + 1.19565217
    
    UH1(ord_uh1, x4 * 24.0 / tstep, D)
    UH2(ord_uh2, x4 * 24.0 / tstep, D)
    
    p_in = 0.0
    epot = 0.0
    q_out = 0.0
    exch = 0.0
    
    Gr4j(st, st_uh1, st_uh2, ord_uh1, ord_uh2, x1, x2, x3, x4, 
    p_in, epot, q_out, exch, tstep, time)
    
end


function get_param_ranges(model::Gr4j)
    
    param_range = Dict(:x1 => (1.0, 20000.0),
    :x2 => (-100.0, 100.0),
    :x3 => (1.0, 20000.0),
    :x4 => (0.5, 10.0))
    
end


function init_states!(model::Gr4j, init_time::DateTime)
    
    model.time = init_time
    
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


function get_water_stored(model::Gr4j)

    water_stored = sum(model.st) + sum(model.st_uh1) + sum(model.st_uh2)

    return water_stored

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
        ER = St[1] * (2.0 - Sr) * TWS / (1.0 + (1.0 - Sr) * TWS)
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
        PS = A * (1.0 - Sr * Sr) * TWS / (1.0 + Sr * TWS)
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
    PERC = St[1] * (1.0 - 1.0 / sqrt(sqrt(1.0 + Sr / scale_param)))
    
    St[1] = St[1] - PERC
    
    PR = PR + PERC
    
    # Split of effective rainfall into the two routing components
    
    PRHU1 = PR * B
    PRHU2 = PR * (1.0 - B)
    
    # Convolution of unit hydrograph UH1
    
    NH = length(OrdUH1)
    
    for K = 1:NH - 1
        StUH1[K] = StUH1[K + 1] + OrdUH1[K] * PRHU1
    end
    
    StUH1[NH] = OrdUH1[NH] * PRHU1
    
    # Convolution of unit hydrograph UH2
    
    NH = length(OrdUH2)
    
    for K = 1:NH - 1
        StUH2[K] = StUH2[K + 1] + OrdUH2[K] * PRHU2
    end
    
    StUH2[NH] = OrdUH2[NH] * PRHU2
    
    # Potential intercatchment semi-exchange
    
    scale_param = (tstep / 24.0)
    
    Rr = St[2] / model.x3
    EXCH = scale_param * model.x2 * Rr * Rr * Rr * sqrt(Rr)
    
    # Routing store
    
    AEXCH1 = EXCH
    
    if (St[2] + StUH1[1] + EXCH) < 0.0
        AEXCH1 = -St[2] - StUH1[1]
    end
    
    St[2] = St[2] + StUH1[1] + EXCH
    
    if St[2] < 0.0
        St[2] = 0.0
    end
    
    scale_param = (24.0 / tstep)
    Rr = St[2]^4 / (model.x3^4 * scale_param)
    
    QR = St[2] * (1.0 - 1.0 / sqrt(sqrt(1.0 + Rr)))
    
    St[2] = St[2] - QR
    
    # Runoff from direct branch QD
    
    AEXCH2 = EXCH
    
    if (StUH2[1] + EXCH) < 0.0
        AEXCH2 = -StUH2[1]
    end
    
    QD = max(0.0, StUH2[1] + EXCH)
    
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
    model.exch    = AEXCH1 + AEXCH2
    
    return nothing
    
end


function SS1(I, C, D)
    
    FI = I
    if FI <= 0.0
        SS1 = 0.0
    elseif FI < C
        SS1 = (FI / C)^D
    else
        SS1 = 1.0
    end
    
end


function SS2(I, C, D)
    
    FI = I
    if FI <= 0.0
        SS2 = 0.0
    elseif FI <= C
        SS2 = 0.5 * (FI / C)^D
    elseif FI < 2.0 * C
        SS2 = 1.0 - 0.5 * (2.0 - FI / C)^D
    else
        SS2 = 1.0
    end
    
end


function UH1(OrdUH1, C, D)
    
    NH = length(OrdUH1)
    
    for I in 1:NH
        
        OrdUH1[I] = SS1(I, C, D) - SS1(I - 1, C, D)
        
    end
    
end


function UH2(OrdUH2, C, D)
    
    NH = length(OrdUH2)
    
    for I in 1:NH
        
        OrdUH2[I] = SS2(I, C, D) - SS2(I - 1, C, D)
        
    end
    
end