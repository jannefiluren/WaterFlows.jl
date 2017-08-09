# Glacier component

mutable struct NoGlacier <: AbstractGlacier
    
    tair::Array{Float64,1}
    q_out::Array{Float64,1}
    frac::Array{Float64,1}
    tstep::Float64
    time::DateTime
    
end


function NoGlacier(tstep::Float64, time::DateTime, frac::Array{Float64,1})
    
    tair  = zeros(Float64, length(frac))
    q_out = zeros(Float64, length(frac))
    
    NoGlacier(tair, q_out, frac, tstep, time)
    
end


function get_param_ranges(model::NoGlacier)
    
    param_range = Dict()
    
end


function init_states!(model::NoGlacier)

    return nothing

end


function run_timestep(m::NoGlacier, s::AbstractSnow)
    
    for reg in eachindex(m.frac)
        
        m.q_out[reg] = 0.0
        
    end
    
    m.time += Dates.Hour(m.tstep)
    
    return nothing
    
end
