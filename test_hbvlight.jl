
using Vann2

# Input data

path = joinpath(Pkg.dir("Vann2"), "data", "atnasjo")

date, tair, prec, q_obs, frac = load_data(path)

date_start = DateTime(2010, 01, 01)

date_stop = DateTime(2015, 01, 01)

date, tair, prec, q_obs = crop_data(date, tair, prec, q_obs, date_start, date_stop)

epot = epot_zero(date)

input = InputPTE(prec, tair, epot)

# Setup model

tstep = 24.0

time = date[1]

lake = 0.0

area = reshape(frac, 1, length(frac))

hydro = HbvLight(tstep, time, area, lake)

model = FullModel(hydro)

# Run model

q_sim = run_model(model, input)

# Run calibration

param_init = get_params(model)

run_model_calib(model, input, q_sim, warmup = 1, verbose = :verbose)

