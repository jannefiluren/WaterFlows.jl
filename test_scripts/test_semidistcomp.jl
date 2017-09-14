
# Load packages

using Vann2
using DataFrames
using PyPlot


# Input data

path = joinpath(Pkg.dir("Vann2"), "data", "fetvatn")

date, tair, prec, q_obs, frac_lus, frac_area, elev = loaddata(path)

date_start = DateTime(2010, 01, 01)

date_stop = DateTime(2015, 01, 01)

date, tair, prec, q_obs = crop_data(date, tair, prec, q_obs, date_start, date_stop)

lat = 70.0

epot = oudin(date, tair, lat, frac_area)

input = InputPTE(prec, tair, epot)


# Test SemiDistComp

tstep = 24.0

time = date[1]

lake = 0.0

snow = HbvLightSnow(tstep, time, frac_lus)

glacier = NoGlacier()

hydro = Hbv(tstep, time)

model = SemiDistComp(snow, glacier, hydro)

q_obs = run_model(model, input)

param_init = get_params(model)

param_tuned = run_model_calib(model, input, q_obs, warmup = 1, verbose = :verbose)

println(round.(param_init,1))

println(round.(param_tuned,1))

set_params!(model, param_tuned)

q_sim = run_model(model, input)

pygui(true)

plot(date, q_sim)
plot(date, q_obs)



