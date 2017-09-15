

using VannModels
using JLD

# Input data

path = joinpath(Pkg.dir("VannModels"), "data", "fetvatn")

date, tair, prec, q_obs, frac_lus, frac_area, elev = load_data(path)

date_start = DateTime(2010, 01, 01)

date_stop = DateTime(2015, 01, 01)

date, tair, prec, q_obs = crop_data(date, tair, prec, q_obs, date_start, date_stop)















