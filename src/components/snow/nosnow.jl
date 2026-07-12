# Snow component

mutable struct NoSnow <: AbstractSnow

    swe::Array{Float64, 2}
    p_in::Array{Float64, 1}
    tair::Array{Float64, 1}
    q_out::Array{Float64, 2}
    frac_lus::Array{Float64, 2}
    tstep::Float64
    time::DateTime

end


function NoSnow(tstep::Float64, time::DateTime, frac_lus::DataFrame)

    @assert (1.0 <= tstep <= 24.0) && isinteger(tstep) "Time step must be a whole number of hours in range 1.0 - 24.0"

    frac_lus = Matrix{Float64}(frac_lus)
    frac_lus = transpose(frac_lus)

    swe = zeros(Float64, size(frac_lus))
    p_in = zeros(Float64, size(frac_lus, 2))
    tair = zeros(Float64, size(frac_lus, 2))
    q_out = zeros(Float64, size(frac_lus))

    return NoSnow(swe, p_in, tair, q_out, frac_lus, tstep, time)

end


function get_param_ranges(model::NoSnow)

    return param_range = Dict()

end


function init_states!(model::NoSnow, init_time::DateTime)

    model.time = init_time

    return nothing

end


function get_water_stored(model::NoSnow)

    return 0.0

end


function run_timestep(m::NoSnow)

    for ireg in 1:size(m.frac_lus, 2)

        for ilus in 1:size(m.frac_lus, 1)

            m.q_out[ilus, ireg] = m.p_in[ireg]

        end

    end

    m.time += Dates.Hour(m.tstep)

    return nothing

end
