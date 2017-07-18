"""
    epot_monthly(datevec)

Use montly average values for potential evapotranspiration.
See report "The nordic HBV model" by Nils Roar Saelthun.

"""

function epot_monthly(datevec)

    # Monthly values of potential evapotranspiration (mm/day)

    epot_month = [0.7, 0.7, 0.7, 1.0, 1.3, 1.4, 1.3, 1.1, 1.0, 0.9, 0.7, 0.7]

    # Assign montly values to days

    epot = zeros(Float64, length(datevec))

    for i in eachindex(datevec)

        imonth = Dates.month(datevec[i])

        epot[i] = epot_month[imonth]

    end

    return(epot)

end


"""
    epot_zero(datevec)

Set potential evapotranspiration to zero.

"""

function epot_zero(datevec)

    epot = zeros(Float64, length(datevec))

    return(epot)

end
