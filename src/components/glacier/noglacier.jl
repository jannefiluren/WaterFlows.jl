# Glacier component

mutable struct NoGlacier <: AbstractGlacier
    
    tair::Array{Float64,1}
    q_out::Array{Float64,1}
    frac::Array{Float64,1}
    iglacier::Int64
    tstep::Float64
    time::DateTime
    
end


function NoGlacier(tstep::Float64, time::DateTime, frac_lus::DataFrame)

    iglacier = find(names(frac_lus) .== :glacier)
    iglacier = iglacier[1]
    
    frac = convert(Array{Float64,2}, frac_lus)
    frac = frac[:, iglacier]
    
    tair  = zeros(Float64, length(frac))
    q_out = zeros(Float64, length(frac))
    
    NoGlacier(tair, q_out, frac, iglacier, tstep, time)
    
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
