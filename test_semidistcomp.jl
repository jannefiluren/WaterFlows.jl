
# Load packages

using Vann2


# Input data

path = joinpath(Pkg.dir("Vann2"), "data", "atnasjo")

date, tair, prec, q_obs, frac = load_data(path)

date_start = DateTime(2010, 01, 01)

date_stop = DateTime(2015, 01, 01)

date, tair, prec, q_obs = crop_data(date, tair, prec, q_obs, date_start, date_stop)

epot = epot_zero(date)

input = InputPTE(prec, tair, epot)


# Test SemiDistComp

tstep = 24.0

time = date[1]

snow = TinBasic(tstep, time, frac)

hydro = Hbv(tstep, time)

model = SemiDistComp(snow, hydro)

q_sim = run_model(model, input)

param_init = get_params(model)

param_out = run_model_calib(model, input, q_sim, warmup = 1, verbose = :verbose)

println(param_init)

println(param_out)


# Test SemiDistFull

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

println(param_out)
