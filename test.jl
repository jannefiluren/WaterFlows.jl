

using VannModels

path = joinpath(Pkg.dir("VannModels"), "data", "atnasjo")

date, tair, prec, q_obs, frac_lus, frac_area, elev = load_data(path)

lat = 60.0

epot = oudin(date, tair, lat, frac_area)

input = InputPTE(prec, tair, epot)

tstep = 24.0

time = date[1]

model = model_hbv_light(tstep, time, frac_lus)

q_sim = run_model(model, input)


snow = HbvLightSnow(tstep, time, frac_lus)

glacier = NoGlacier()

subsurf = Gr4j(tstep, time)

model = ModelComp(snow, glacier, subsurf)

q_sim = run_model(model, input)
