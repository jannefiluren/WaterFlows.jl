using VannModels

path = joinpath(@__DIR__, "data", "atnasjo")
date, tair, prec, q_obs, frac_lus, frac_area, elev = load_data(path)
lat = 60.0
epot = oudin(date, tair, lat, frac_area)
input = InputPTE(date, prec, tair, epot)
tstep = 24.0
tstart = date[1]
model = model_hbv_light(tstep, tstart, frac_lus)
q_sim = run_model(model, input)

snow = HbvLightSnow(tstep, tstart, frac_lus)
glacier = NoGlacier()
subsurf = Gr4j(tstep, tstart)
model = ModelComp(snow, glacier, subsurf)
q_sim = run_model(model, input)
param_tuned = run_model_calib(model, input, q_obs, warmup = 1, verbose = :silent)
set_params!(model, param_tuned)
q_sim = run_model(model, input)
