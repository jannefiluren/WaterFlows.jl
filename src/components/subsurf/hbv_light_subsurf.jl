# Hydrological component
#
# Rate parameters are given in daily units and are rescaled
# to the model time step internally in run_timestep.

mutable struct HbvLightSubsurf <: AbstractSubsurfDist

    sm::Array{Float64, 2}        # Soil storage [mm]
    suz::Float64                 # Storage upper zone [mm]
    slz::Float64                 # Storage lower zone [mm]
    st_uh::Array{Float64, 1}     # Storage unit hydrograph [mm]
    ord_uh::Array{Float64, 1}    # Ordinates of unit hydrograph [-]
    perc::Float64                # Maximum flow from upper to lower zone [mm d^-1]
    k0::Float64                  # Recession coefficient upper zone above uzl threshold [d^-1]
    k1::Float64                  # Recession coefficient upper zone [d^-1]
    k2::Float64                  # Recession coefficient lower zone [d^-1]
    uzl::Float64                 # Threshold for fast runoff from upper zone [mm]
    fc::Array{Float64, 1}        # Maximum of soil storage [mm]
    lp::Array{Float64, 1}        # Fraction of fc below which evaporation is reduced [-]
    beta::Array{Float64, 1}      # Shape coefficient [-]
    maxbas::Float64              # Routing, length of weighting function [d]
    snow::Array{Bool, 2}         # Mask indicating presence of snow [-]
    p_in::Array{Float64, 2}      # Input from precipitation and snowmelt [mm d^-1]
    epot::Float64                # Potential evapotranspiration [mm d^-1]
    q_out::Float64               # Runoff [mm d^-1]
    aevap::Float64               # Actual evapotranspiration [mm d^-1]
    frac_lus::Array{Float64, 2}  # Landuse fractions [-]
    tstep::Float64               # Time step [h]
    time::DateTime               # Current time [-]

end


function HbvLightSubsurf(tstep::Float64, time::DateTime, frac_lus::DataFrame)

    @assert (1.0 <= tstep <= 24.0) && isinteger(tstep) "Time step must be a whole number of hours in range 1.0 - 24.0"

    frac_lus = Matrix{Float64}(frac_lus)
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
    aevap = 0.0

    ord_uh = compute_hbv_ord(maxbas, tstep)
    st_uh = zero(ord_uh)

    return HbvLightSubsurf(
        sm, suz, slz, st_uh, ord_uh, perc, k0, k1, k2, uzl, fc,
        lp, beta, maxbas, snow, p_in, epot, q_out, aevap, frac_lus, tstep, time
    )

end


function get_param_ranges(m::HbvLightSubsurf)

    return param_range = Dict(
        :perc => (0.1, 1000.0),
        :k0 => (0.001, 0.999),
        :k1 => (0.001, 0.999),
        :k2 => (0.001, 0.999),
        :uzl => (1.0, 1000.0),
        :fc => (50.0, 500.0),
        :lp => (0.3, 1.0),
        :beta => (1.0, 6.0),
        :maxbas => (1.0, 7.0)
    )

end


function init_states!(m::HbvLightSubsurf, init_time::DateTime)

    m.time = init_time

    m.sm .= zero(m.sm)
    m.suz = 0.0
    m.slz = 0.0

    m.ord_uh = compute_hbv_ord(m.maxbas, m.tstep)

    m.st_uh = zero(m.ord_uh)

    return nothing

end


function get_water_stored(m::HbvLightSubsurf)

    water_stored = sum(m.sm .* m.frac_lus) + m.suz + m.slz + sum(m.st_uh)

    return water_stored

end


function run_timestep(m::HbvLightSubsurf)

    # Scale parameters from daily values to the current time step
    dt = m.tstep / 24.0
    k0 = 1.0 - (1.0 - m.k0)^dt
    k1 = 1.0 - (1.0 - m.k1)^dt
    k2 = 1.0 - (1.0 - m.k2)^dt
    perc = m.perc * dt

    epot = m.epot

    sum_to_gw = 0.0
    avg_aet = 0.0

    # Soil routine
    for ireg in 1:size(m.frac_lus, 2)

        for ilus in 1:size(m.frac_lus, 1)

            if m.frac_lus[ilus, ireg] > 0.0

                sm = m.sm[ilus, ireg]
                fc = m.fc[ilus]
                lp = m.lp[ilus]
                beta = m.beta[ilus]

                insoil = m.p_in[ilus, ireg]

                to_gw = 0.0
                old_sm = sm

                # Compute infiltration to soil and recharge to groundwater
                if insoil > 0.0
                    if insoil < 1.0
                        insoil_resid = insoil
                    else
                        # Compute fractions in steps to account for non-linearities
                        insoil_floored = floor(insoil)
                        insoil_resid = insoil - insoil_floored
                        for _ in 1:insoil_floored
                            frac_to_gw = (sm / fc)^beta
                            if frac_to_gw > 1.0
                                frac_to_gw = 1.0
                            end
                            sm = sm + 1.0 - frac_to_gw
                            to_gw = to_gw + frac_to_gw
                        end
                    end
                    frac_to_gw = (sm / fc)^beta
                    if frac_to_gw > 1.0
                        frac_to_gw = 1.0
                    end
                    sm = sm + (1 - frac_to_gw) * insoil_resid
                    to_gw = to_gw + frac_to_gw * insoil_resid
                end

                mean_sm = (sm + old_sm) / 2.0
                if mean_sm < (lp * fc)
                    aet = epot * mean_sm / (lp * fc)
                else
                    aet = epot
                end
                #if snow        # Currently snow does not influence actual evapotranspiration
                #    aet = 0.0
                #end
                sm = sm - aet
                if sm < 0.0
                    sm = 0.0
                end

                avg_aet = avg_aet + aet * m.frac_lus[ilus, ireg]
                sum_to_gw = sum_to_gw + to_gw * m.frac_lus[ilus, ireg]
                m.sm[ilus, ireg] = sm

            end

        end

    end

    # Add groundwater recharge to upper zone storage
    m.suz = m.suz + sum_to_gw

    # Handle flow from upper to lower storage zone
    if (m.suz - perc) < 0.0
        # Move all upper zone water to lower zone if less than percolation rate
        m.slz = m.slz + m.suz
        m.suz = 0.0
    else
        # Remove percolation from upper and add to lower zone storage
        m.slz = m.slz + perc
        m.suz = m.suz - perc
    end

    # Compute outflow from upper zone
    q_suz1 = k1 * m.suz
    if m.suz < m.uzl
        q_suz0 = 0.0
    else
        q_suz0 = k0 * (m.suz - m.uzl)
    end

    # Limit outflow to available water in upper storage zone
    q_suz = min(q_suz1 + q_suz0, m.suz)

    # Update upper zone storage
    m.suz = m.suz - q_suz

    # Outflow from lower zone
    q_slz = k2 * m.slz

    # Update lower zone storage
    m.slz = m.slz - q_slz

    # Total outflow to unit hydrograph
    q_gen = q_suz + q_slz

    # Convolution of unit hydrograph
    nh = length(m.ord_uh)
    for k in 1:(nh - 1)
        m.st_uh[k] = m.st_uh[k + 1] + m.ord_uh[k] * q_gen
    end
    m.st_uh[nh] = m.ord_uh[nh] * q_gen

    # Output runoff
    m.q_out = m.st_uh[1]
    m.st_uh[1] = 0

    # Output actual evapotranspiration
    m.aevap = avg_aet

    # Update time
    m.time += Dates.Hour(m.tstep)

    return nothing

end
