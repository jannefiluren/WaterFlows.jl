using Distributions

maxbas = 3.0

triang = Distributions.TriangularDist(0, maxbas)

# Unit hydrograph ordinates for a triangular weighting function, sampled at
# two time steps and with two sampling conventions:
#
#   Left edge sampling: cdf(i * dt)         -> ordinate k covers [(k-1)*dt, k*dt]
#   Midpoint sampling:  cdf((i - 0.5) * dt) -> ordinate k centred on delay (k-1)*dt
#
# The convolution releases ordinate k with a delay of (k-1)*dt, so the mean
# delay (centroid) is what determines the timing of the routed hydrograph.


# Ordinates for a given time step and sampling convention (edge = 0.0 or 0.5)

function compute_ord(triang, maxbas, tstep, edge)

    dt = tstep / 24.0

    nord = ceil(Int, maxbas / dt)

    cdf_vals = [Distributions.cdf(triang, (i - edge) * dt) for i in 0:nord]

    ord_uh = diff(cdf_vals)

    ord_uh ./= sum(ord_uh)

    return ord_uh, dt

end


# Mean delay (centroid) of the unit hydrograph, i.e. how the convolution times it

centroid(ord_uh, dt) = sum(ord_uh[k] * (k - 1) * dt for k in eachindex(ord_uh))


# Left edge sampling

daily_ord_uh, daily_dt = compute_ord(triang, maxbas, 24.0, 0.0)

hourly_ord_uh, hourly_dt = compute_ord(triang, maxbas, 1.0, 0.0)

# ... produces same mass over the first day

daily_ord_uh[1]

sum(hourly_ord_uh[1:24])

# ... but the centroids differ by ~half a day (the shift)

centroid(daily_ord_uh, daily_dt)

centroid(hourly_ord_uh, hourly_dt)


# Midpoint sampling

daily_ord_uh, daily_dt = compute_ord(triang, maxbas, 24.0, 0.5)

hourly_ord_uh, hourly_dt = compute_ord(triang, maxbas, 1.0, 0.5)

# ... produces mass over the first day that no longer matches

daily_ord_uh[1]

sum(hourly_ord_uh[1:24])

# ... but the centroids now agree (shift-free)

centroid(daily_ord_uh, daily_dt)

centroid(hourly_ord_uh, hourly_dt)
