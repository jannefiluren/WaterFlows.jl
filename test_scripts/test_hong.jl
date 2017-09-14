
using Vann2
using Base.Test

# Input data

path = joinpath(Pkg.dir("Vann2"), "data", "atnasjo")

date, tair, prec, q_obs, frac = loaddata(path)

input = InputPT(prec, tair)

# Setup model

tstep = 24.0

time = date[1]

snow = TinSnow(tstep, time, frac)

model = SnowDistModel(snow)

# Run model

swe_sim = run_model(model, input)

# Run calibration

param_ref = get_params(model)

run_model_calib(model, input, swe_sim, warmup = 1)




