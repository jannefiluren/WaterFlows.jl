"""
Inverse relative distance between earth to sun from day-of-year.
"""
function inverse_dist(doy)
    dr = 1 + 0.033cos(2π/365.0*doy)
end


"""
Solar declination from day-of-year.
"""
function solar_decl(doy)
    δ = 0.409 * sin(2π/365.0*doy - 1.39)
end

    
"""
Sunset hour angle from latitude (rad) and solar declination (rad).
"""
function sunset_hour_angle(lat, δ)
    tmp = -tan(δ)*tan(lat)
    tmp = tmp > 1.0 ? 1.0 : tmp
    tmp = tmp < -1.0 ? -1.0 : tmp
    ω_s = acos(tmp)
end


"""
Maximum daylight hours (hours) from sunset hour angle (rad).
"""
function daylight_hours(ω_s)
    N = (24.0/π)*ω_s
end


"""
Extraterrestrial radiation (MJ m-2 day-1) from date and latitude (°).
"""
function extra_ter_rad(date, lat)
    lat = deg2rad(lat)
    doy = Dates.dayofyear(date)
    dr = inverse_dist(doy)
    δ = solar_decl(doy)
    ω_s = sunset_hour_angle(lat, δ)
    R_a = 1440/π*0.082*dr*(ω_s*sin(lat)*sin(δ) + cos(lat)*cos(δ)*sin(ω_s))
end


"""
Clear-sky solar radiation (MJ m-2 day-1) from date, latitude (°) and
elevation (m).
"""
function clear_sky_rad(date, lat, elev)
    R_a = extra_ter_rad(date, lat)
    R_so = (0.75 + 2e-5*elev)*R_a    
end