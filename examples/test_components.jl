
# Load packages

using WaterFlows
using DataFrames
using GLMakie
using Dates

# Input data

path = joinpath(dirname(pathof(WaterFlows)), "..", "data", "atnasjo")

date, tair, prec, q_obs, frac_lus, frac_area, elev = load_data(path)

date_start = DateTime(2010, 01, 01)

date_stop = DateTime(2015, 01, 01)

date, tair, prec, q_obs = crop_data(date, tair, prec, q_obs, date_start, date_stop)

lat = 70.0

epot = oudin(date, tair, lat, frac_area)

input = InputPTE(date, prec, tair, epot)


# Test model components

tstep = 24.0

time = date[1]

lake = 0.0

snow = HbvLightSnow(tstep, time, frac_lus)

glacier = NoGlacier()

subsurf = Hbv(tstep, time)

model = ModelComp(snow, glacier, subsurf)

init_states!(model, input.time[1])

q_ref = run_model(model, input)

param_init = get_params(model)

param_tuned = run_model_calib(model, input, q_ref, warmup = 1, verbose = :verbose, max_steps = 50000)

println(round.(param_init, digits=1))

println(round.(param_tuned, digits=1))

set_params!(model, param_tuned)

init_states!(model, input.time[1])

q_sim = run_model(model, input)

f = Figure()
ax = Axis(f[1, 1], ylabel="Runoff [mm d^-1]")
lines!(ax, date, q_ref, label="Reference (synthetic)")
lines!(ax, date, q_sim, label="Simulated", linestyle=:dash)
axislegend(ax)
isinteractive() ? display(f) : wait(display(f))
