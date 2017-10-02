# Glacier component

mutable struct NoGlacier <: AbstractGlacier
    
    #= tair::Array{Float64,1}
    q_out::Array{Float64,1}
    frac_lus::Array{Float64,1}
    iglacier::Int64
    tstep::Float64
    time::DateTime =#
    
end


#= function NoGlacier(tstep::Float64, time::DateTime, frac_lus::DataFrame)

    iglacier = find(names(frac_lus) .== :glacier)
    iglacier = iglacier[1]
    
    frac_lus = convert(Array{Float64,2}, frac_lus)
    frac_lus = frac_lus[:, iglacier]
    
    tair  = zeros(Float64, length(frac_lus))
    q_out = zeros(Float64, length(frac_lus))
    
    NoGlacier(tair, q_out, frac_lus, iglacier, tstep, time)
    
end =#


function get_param_ranges(model::NoGlacier)
    
    param_range = Dict()
    
end


function init_states!(model::NoGlacier, init_time::DateTime)

    return nothing

end


function run_timestep(m::NoGlacier, s::AbstractSnow)

    return nothing

end

#= function run_timestep(m::NoGlacier, s::AbstractSnow)
    
    for reg in eachindex(m.frac_lus)
        
        m.q_out[reg] = 0.0
        
    end
    
    m.time += Dates.Hour(m.tstep)
    
    return nothing
    
end
 =#