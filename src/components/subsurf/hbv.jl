# Hydrological component
#
# Rate parameters are given in daily units and are rescaled
# to the model time step internally in run_timestep.

mutable struct Hbv <: AbstractSubsurfLumped
    
    sm::Float64                 # Soil storage [mm]
    suz::Float64                # Storage upper zone [mm]
    slz::Float64                # Storage lower zone [mm]
    st_uh::Array{Float64,1}     # Storage unit hydrograph [mm]
    ord_uh::Array{Float64,1}    # Ordinates of unit hydrograph [-]
    fc::Float64                 # Maximum of soil storage [mm]
    lp::Float64                 # Fraction of fc below which evaporation is reduced [-]
    k0::Float64                 # Recession coefficient upper zone above ulz threshold [d^-1]
    k1::Float64                 # Recession coefficient upper zone [d^-1]
    k2::Float64                 # Recession coefficient lower zone [d^-1]
    beta::Float64               # Shape coefficient [-]
    perc::Float64               # Maximum flow from upper to lower zone [mm d^-1]
    ulz::Float64                # Threshold for fast runoff from upper zone [mm]
    maxbas::Float64             # Routing, length of weighting function [d]
    p_in::Float64               # Input from precipitation and snowmelt [mm d^-1]
    epot::Float64               # Potential evapotranspiration [mm d^-1]
    q_out::Float64              # Runoff [mm d^-1]
    aevap::Float64              # Actual evapotranspiration [mm d^-1]
    tstep::Float64              # Time step [h]
    time::DateTime              # Current time [-]
    
end


function Hbv(tstep::Float64, time::DateTime)
    
    @assert (1.0 <= tstep <= 24.0) && isinteger(tstep) "Time step must be a whole number of hours in range 1.0 - 24.0"

    sm     = 0.0
    suz    = 0.0
    slz    = 0.0

    fc, lp, k0, k1, k2, beta, perc, ulz, maxbas = (100.0, 0.8, 0.05, 0.05, 0.01, 1.0, 2.0, 30., 2.5)

    p_in = 0.0
    epot = 0.0
    q_out = 0.0
    aevap = 0.0

    ord_uh = compute_hbv_ord(maxbas, tstep)
    st_uh  = zero(ord_uh)
    
    Hbv(sm, suz, slz, st_uh, ord_uh, fc, lp, k0, k1, k2, beta, perc, ulz, maxbas, p_in, epot, q_out, aevap, tstep, time)
    
end


function get_param_ranges(m::Hbv)
    
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


function init_states!(m::Hbv, init_time::DateTime)
    
    m.time = init_time
    
    m.sm      = 0.5*m.fc
    m.suz     = 0.5*m.ulz
    m.slz     = 0.5*m.ulz
    
    m.ord_uh = compute_hbv_ord(m.maxbas, m.tstep)

    m.st_uh = zero(m.ord_uh)

end


function get_water_stored(m::Hbv)

    water_stored = m.sm + m.suz + m.slz + sum(m.st_uh)

    return water_stored

end


function run_timestep(m::Hbv)

    # Scale parameters from daily values to the current time step

    dt = m.tstep / 24.0

    k0 = 1.0 - (1.0 - m.k0)^dt
    k1 = 1.0 - (1.0 - m.k1)^dt
    k2 = 1.0 - (1.0 - m.k2)^dt

    perc = m.perc * dt

    # Soil moisture zone (assume no evaporation during rainfall)
    
    if m.p_in > 0.0
        
        # Beta function
        
        f_recharge = (m.sm / m.fc) ^ m.beta
        
        # Groundwater recharge
        
        recharge = f_recharge * m.p_in
        
        # Update soil moisture zone
        
        m.sm = m.sm + m.p_in - recharge
        
        # Add excess soil moisture to groundwater recharge
        
        if m.sm > m.fc
            recharge += m.sm - m.fc
            m.sm = m.fc
        end
        
        # No evapotranspiration
        
        aevap = 0.0
        
    else
        
        # Compute actual evapotranspiration
        
        aevap = m.epot * min(m.sm/(m.fc*m.lp), 1.0)
        
        # Update soil moisture zone
        
        m.sm = m.sm - aevap
        
        # Check limits for soil moisture zone
        
        if m.sm < 0.0
            aevap = max(aevap + m.sm, 0.0)
            m.sm = 0.0
        end
        
        # No groundwater recharge
        
        recharge = 0.0
        
    end
    
    # Add recharge to upper groundwater box
    
    m.suz = m.suz + recharge
    
    # Remove percolation from upper groundwater box
    
    perc_now = min(perc, m.suz)
    
    m.suz = m.suz - perc_now
    
    # Compute runoff from upper groundwater box and update storage
    
    q_suz = min(k1 * m.suz + k0 * max(m.suz-m.ulz, 0.0), m.suz)

    m.suz = m.suz - q_suz
    
    # Add precolation to lower groundwater box
    
    m.slz = m.slz + perc_now
    
    # Compute runoff from lower groundwater box and update storage
    
    q_slz = k2 * m.slz
    
    m.slz = m.slz - q_slz
    
    # Convolution of unit hydrograph
    
    q_tmp = q_suz + q_slz
    
    nh = length(m.ord_uh)
    
    for k = 1:nh-1
        m.st_uh[k] = m.st_uh[k+1] + m.ord_uh[k]*q_tmp
    end
    
    m.st_uh[nh] = m.ord_uh[nh] * q_tmp
    
    # Compute total runoff
    
    m.q_out = m.st_uh[1]
    m.st_uh[1] = 0

    # Add actual evapotranspiration

    m.aevap = aevap

    # Update time
    
    m.time  += Dates.Hour(m.tstep)
    
    return nothing
    
end


function compute_hbv_ord(maxbas, tstep = 24.0)

    # Unit hydrograph ordinates from a triangular weighting function
    # with base length maxbas [d] sampled at the model time step

    dt = tstep / 24.0

    nord = ceil(Int, maxbas / dt)

    triang = Distributions.TriangularDist(0, maxbas)
    triang_cdf = [Distributions.cdf(triang, i * dt) for i in 0:nord]
    ord_uh = diff(triang_cdf)

end

