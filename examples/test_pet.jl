# Load packages

using WaterFlows
using GLMakie

# Input data

path = joinpath(dirname(pathof(WaterFlows)), "..", "data", "atnasjo")

date, tair, prec, q_obs, frac_lus, frac_area = load_data(path)

lat = 48.0   # Latitude in degrees (oudin and hamon convert to radians internally)

pet_oudin = oudin(date, tair, lat, frac_area)

pet_hamon = hamon(date, tair, lat, frac_area)

f = Figure()
ax = Axis(f[1, 1], ylabel = "Potential evapotranspiration [mm d^-1]")
lines!(ax, date, pet_oudin, label = "Oudin")
lines!(ax, date, pet_hamon, label = "Hamon")
axislegend(ax)
isinteractive() ? display(f) : wait(display(f))
