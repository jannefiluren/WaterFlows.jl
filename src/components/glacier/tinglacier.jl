# Glacier component

mutable struct TinGlacier <: AbstractGlacier

    ddf::Float64
    tair::Array{Float64,1}
    q_out::Array{Float64,1}
    frac_lus::Array{Float64,1}
    iglacier::Int64
    tstep::Float64
    time::DateTime
    
end


function TinGlacier(tstep::Float64, time::DateTime, frac_lus::DataFrame)

    @assert (1.0 <= tstep <= 24.0) && isinteger(tstep) "Time step must be a whole number of hours in range 1.0 - 24.0"

    iglacier = findfirst(==("glacier"), names(frac_lus))

    @assert iglacier !== nothing "frac_lus must contain a column named glacier"

    frac_lus = Matrix{Float64}(frac_lus)
    frac_lus = frac_lus[:, iglacier]
    
    tair  = zeros(Float64, length(frac_lus))
    q_out = zeros(Float64, length(frac_lus))

    ddf = 5.0
    
    TinGlacier(ddf, tair, q_out, frac_lus, iglacier, tstep, time)
    
end


function get_param_ranges(model::TinGlacier)
    
    param_range = Dict(:ddf => (0.1, 10.0))
    
end


function init_states!(model::TinGlacier, init_time::DateTime)
    
    model.time = init_time

    return nothing

end


function run_timestep(g::TinGlacier, s::AbstractSnow)

    # Scale degree-day factor from daily value to the current time step

    dt = g.tstep / 24.0

    for reg in eachindex(g.frac_lus)

        if s.swe[g.iglacier, reg] > 0.0
            g.q_out[reg] = 0.0
        else
            g.q_out[reg] = dt * pot_melt(g.tair[reg], g.ddf)
        end
                
    end
    
    g.time += Dates.Hour(g.tstep)
    
    return nothing
    
end
