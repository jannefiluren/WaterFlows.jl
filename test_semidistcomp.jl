
# Load packages

using Vann2
using DataFrames

# Input data

path = joinpath(Pkg.dir("Vann2"), "data", "fetvatn")

date, tair, prec, q_obs, frac_lus, frac_area, elev = load_data(path)

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

glacier = HockGlacier(tstep, time, frac_lus, lat, elev)

hydro = HbvLightSubsurf(tstep, time, frac_lus, lake)

model = SemiDistComp(snow, glacier, hydro)

q_sim = run_model(model, input)

param_init = get_params(model)

param_out = run_model_calib(model, input, q_sim, warmup = 1, verbose = :verbose)

println(round.(param_init,2))

println(round.(param_out,2))






#=
# Test SemiDistComp

tstep = 24.0

time = date[1]

snow = TinSnow(tstep, time, frac_lus)

glacier = TinGlacier(tstep, time, frac_lus)

hydro = Gr4j(tstep, time)

model = SemiDistComp(snow, glacier, hydro)

q_sim = run_model(model, input)

param_init = get_params(model)

param_out = run_model_calib(model, input, q_sim, warmup = 1, verbose = :verbose)

println(round.(param_init,2))

println(round.(param_out,2))
=#




























#= Test SemiDistFull

tstep = 24.0

time = date[1]

lake = 0.0

area = reshape(frac, 1, length(frac))

hydro = HbvLight(tstep, time, area, lake)

model = SemiDistFull(hydro)

q_sim = run_model(model, input)

param_init = get_params(model)

param_out = run_model_calib(model, input, q_sim, warmup = 1, verbose = :verbose)

println(param_init)

println(param_out)=#
