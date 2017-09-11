# Snow component

mutable struct NoSnow <: AbstractSnow
    
    p_in::Array{Float64,1}
    tair::Array{Float64,1}
    q_out::Array{Float64,1}   # This needs to be a 2d array
    frac::Array{Float64,1}    # This needs to be a 2d array
    tstep::Float64
    time::DateTime
    
end


function NoSnow(tstep::Float64, time::DateTime, frac::Array{Float64,1})
    
    p_in  = zeros(Float64, length(frac))
    tair  = zeros(Float64, length(frac))
    q_out = zeros(Float64, length(frac))
    
    NoSnow(p_in, tair, q_out, frac, tstep, time)
    
end


function get_param_ranges(model::NoSnow)
    
    param_range = Dict()
    
end


function init_states!(model::NoSnow)

    return nothing

end


function run_timestep(m::NoSnow)
    
    for reg in eachindex(m.frac)
        
        m.q_out[reg] = m.p_in[reg]
        
    end
    
    m.time += Dates.Hour(m.tstep)
    
    return nothing
    
end
