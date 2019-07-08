

using WaterFlows
using JLD

############################################################################################

# Load all input data, and select one time period

path = joinpath(Pkg.dir("WaterFlows"), "data", "fetvatn")

date, tair, prec, q_obs, frac_lus, frac_area, elev = load_data(path)

date_start = DateTime(2010, 01, 01)

date_stop = DateTime(2015, 01, 01)

date, tair, prec, q_obs = crop_data(date, tair, prec, q_obs, date_start, date_stop)

lat = 60.0

epot = oudin(date, tair, lat, frac_area)

input = InputPTE(date, prec, tair, epot)

############################################################################################

# Run the model for this period

tstep = 24.0

time = date[1]

model = setup_hbv_light(tstep, time, frac_lus)

q_sim = run_model(model, input)

############################################################################################

# Save the final state of the model to a file

save("model_state.jld", "model", model)

############################################################################################

# Load the states of the model again

tmp = load("model_state.jld")

tmp["model"]

############################################################################################

# Now we could run the model from here again ...

