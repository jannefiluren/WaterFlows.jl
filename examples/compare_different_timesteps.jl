# Load packages

using WaterFlows
using DataFrames
using GLMakie
using Dates

# Construct a subsurface component (constructors differ in their arguments)

setup_subsurf(::Type{Gr4j},            tstep, time, frac_lus) = Gr4j(tstep, time)
setup_subsurf(::Type{Hbv},             tstep, time, frac_lus) = Hbv(tstep, time)
setup_subsurf(::Type{HbvLightSubsurf}, tstep, time, frac_lus) = HbvLightSubsurf(tstep, time, frac_lus)

# Helper function

function run_example(input_daily, input_hourly, frac_lus, snow_model, subsurf_model)

    # Run daily model

    snow = snow_model(24.0, input_daily.time[1], frac_lus)

    glacier = NoGlacier()

    subsurf = setup_subsurf(subsurf_model, 24.0, input_daily.time[1], frac_lus)

    model = ModelComp(snow, glacier, subsurf)

    init_states!(model, input_daily.time[1])

    q_daily = run_model(model, input_daily)

    # Run hourly model

    snow = snow_model(1.0, input_hourly.time[1], frac_lus)

    glacier = NoGlacier()

    subsurf = setup_subsurf(subsurf_model, 1.0, input_hourly.time[1], frac_lus)

    model = ModelComp(snow, glacier, subsurf)

    init_states!(model, input_hourly.time[1])

    q_hourly = run_model(model, input_hourly)

    # Plot results

    q_hourly_daily = [sum(@view q_hourly[(i-1)*24+1 : i*24]) for i in eachindex(input_daily.time)]

    fig = Figure()
    ax = Axis(fig[1, 1], title = "$(snow_model) + $(subsurf_model)", ylabel = "Runoff [mm/day]")
    lines!(ax, time_daily, q_daily, label = "Daily")
    lines!(ax, time_daily, q_hourly_daily, label = "Hourly (aggregated)")
    fig[1, 2] = Legend(fig, ax, "Timestep", framevisible = false)
    display(GLMakie.Screen(), fig)

    return nothing

end

# Daily input data

path = joinpath(dirname(pathof(WaterFlows)), "..", "data", "atnasjo")

time_daily, tair_daily, prec_daily, q_obs, frac_lus, frac_area, elev = load_data(path)

epot_daily = oudin(time_daily, tair_daily, 70.0, frac_area)

input_daily = InputPTE(time_daily, prec_daily, tair_daily, epot_daily)

# Hourly input

time_hourly = [time_daily[j] + Hour(i) for i in 0:23, j in eachindex(time_daily)][:]

time_hourly = DateTime(first(time_daily)):Hour(1):DateTime(last(time_daily)) + Hour(23)
tair_hourly = repeat(tair_daily, inner=(1, 24))
prec_hourly = repeat(prec_daily, inner=(1, 24)) / 24

epot_hourly = repeat(epot_daily, inner=24) / 24

input_hourly = InputPTE(time_hourly, prec_hourly, tair_hourly, epot_hourly)

# Run model
#
# The unit hydrographs are sampled at interval midpoints, so the routing is
# discretised consistently across time steps. The Hbv-based subsurface models
# therefore line up closely between daily and hourly runs. Gr4j still shows a
# residual daily/hourly difference by design: its unit-hydrograph shape exponent
# varies with the time step and its production/routing stores are nonlinear, so
# it is time-step sensitive beyond the routing discretisation.

snow_models = [HbvLightSnow, TinSnow]
subsurf_models = [Gr4j, HbvLightSubsurf, Hbv]

for snow_model in snow_models, subsurf_model in subsurf_models
    run_example(input_daily, input_hourly, frac_lus, snow_model, subsurf_model)
end
