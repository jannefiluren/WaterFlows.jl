


"""Potential evapotranspiration using Oudins formula (mm/day)."""
function oudin(date::Array{DateTime,1}, tair::Array{Float64,2}, lat::Float64, frac_area::Array{Float64,1})

    pet = similar(tair)

    for reg in eachindex(frac_area)
        pet[reg,:] = oudin(date, tair[reg,:] , lat) * frac_area[reg]
    end

    pet = sum(pet, 1)

    squeeze(pet,1)

end


"""Potential evapotranspiration using Oudins formula (mm day-1)."""
function oudin(date::Array{DateTime,1}, tair::Array{Float64,1}, lat::Float64)

    oudin_tmp(tair, date) = oudin(tair, date, lat)

    pet = map(oudin_tmp, tair, date)

end


"""Potential evapotranspiration using Oudins formula (mm day-1)."""
function oudin(tair::Float64, date::DateTime, lat::Float64)

    # Convert date to day-of-year

    doy = Dates.dayofyear(date)

    # Converts latitude in radians

    fi = deg2rad(lat)
    cosfi = cos(fi)
    afi = abs(lat/42.0)

    # Declination of the sun in radians

    teta = 0.4093*sin(doy/58.1-1.405)
    costeta = cos(teta)
    cosgz = max(0.001,cos(fi-teta))

    # Noon angular zenith distance of the sun

    cosgz2 = cosgz*cosgz

    if cosgz2 >= 1.0
        singz = 0.0
    else
        singz = sqrt(1.-cosgz2)
    end

    cosom = 1.0-cosgz/cosfi/costeta

    if cosom < -1.0
        cosom = -1.0
    end

    if cosom > 1.0
        cosom = 1.0
    end

    cosom2 = cosom*cosom

    if cosom2 >= 1.0
        sinom = 0.0
    else
        sinom = sqrt(1.0-cosom2)
    end

    om = acos(cosom)

    # Average angular zenith distance of the sun

    cospz = cosgz+cosfi*costeta*(sinom/om-1.)

    if cospz < 0.001
        cospz = 0.001
    end

    # Radius vector of the sun

    eta = 1.0+cos(doy/58.1)/30.0

    # Extra-atmospheric global radiation 

    ge = 446.0*om*cospz*eta

    # Potential evapotranspiration

    pet = max(0.0, ge*(tair+5.0)/100.0/28.5)

end


"""Potential evapotranspiration using Hamons formula (mm/day)."""
function hamon(date::Array{DateTime,1}, tair::Array{Float64,2}, lat::Float64, frac_area::Array{Float64,1})

    pet = similar(tair)

    for reg in eachindex(frac_area)
        pet[reg,:] = hamon(date, tair[reg,:] , lat) * frac_area[reg]
    end

    pet = sum(pet, 1)

    squeeze(pet,1)

end


"""Potential evapotranspiration using Hamons formula (mm/day)."""
function hamon(date::Array{DateTime,1}, tair::Array{Float64,1}, lat::Float64)

    # Compute monthly average temperature

    tair = monthly_mean(date, tair)

    # Compute monthly average day length

    day_length_tmp(date) = day_length(date, lat)

    ld = map(day_length_tmp, date)

    ld = monthly_mean(date, ld)

    # Compute saturation vapor density

    esat = map(sat_vap_pressure, tair)

    # Compute potential evapotranspiration

    pet = map(hamon, tair, esat, ld)

end


"""Potential evapotranspiration using Hamons formula (mm/day)."""
function hamon(tair::Float64, esat::Float64, ld::Float64)

    pet = 29.8 * ld * esat / (tair + 273.2)

end


"""Compute monthly averages."""
function monthly_mean(date, values)

    months = map(Dates.month, date)
    
    for month in unique(months)
        
        imonth = find(month .== months)
        
        values[imonth] .= mean(values[imonth]) 
        
    end

    return values   
    
end


"""Saturated vapor pressure (kPa)"""
function sat_vap_pressure(tair)
    
    esat = 0.6108*exp(17.27*tair / (tair+237.3))
    
end


"""Solar declination."""
function solar_decl(doy)
    
    0.409 * sin(2π/365.0*doy - 1.39)

end


"""Sunset hour angle."""
function sunset_hour_angle(lat, δ)

    tmp = -tan(δ)*tan(lat)

    tmp = tmp > 1.0 ? 1.0 : tmp
    tmp = tmp < -1.0 ? -1.0 : tmp
    
    acos(tmp)

end


"""Day length (hours)"""
function day_length(date, lat)
    
    # Compute day of year
    
    doy = Dates.dayofyear(date)
    
    # Declination
    
    δ = solar_decl(doy)
    
    # Sunset hour angle
    
    ω = sunset_hour_angle(lat, δ)
    
    # Day length
    
    (24.0/π)*ω
    
end


"""Extraterrestrial radiation (MJ m-2 day-1)"""
function extra_ter_rad(date, lat)

    # Compute day of year
    
    doy = Dates.dayofyear(date)

    # Declination
    
    δ = solar_decl(doy)
    
    # Sunset hour angle
    
    ω = sunset_hour_angle(lat, δ)

    # Inverse relative distance earth to sun

    dr = 1 + 0.033cos(2π/365.0*doy)

    # Extraterrestrial radiation

    1440/π*0.082*dr*(ω*sin(lat)*sin(δ) + cos(lat)*cos(δ)*sin(ω))

end


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
