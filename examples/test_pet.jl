
# Load packages

using WaterFlows
using PyPlot

# Input data

path = joinpath(Pkg.dir("WaterFlows"), "data", "atnasjo")

date, tair, prec, q_obs, frac_lus, frac_area = load_data(path)

lat = 48   # Latitude

lat = deg2rad(lat) 

pet_oudin = oudin(date, tair, lat, frac_area)

pet_hamon = hamon(date, tair, lat, frac_area)

pygui(true)
plot(date, pet_oudin)
plot(date, pet_hamon)
