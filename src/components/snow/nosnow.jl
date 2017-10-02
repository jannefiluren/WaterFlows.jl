# Snow component

mutable struct NoSnow <: AbstractSnow
    
    p_in::Array{Float64,1}
    tair::Array{Float64,1}
    q_out::Array{Float64,1}   # This needs to be a 2d array
    frac_lus::Array{Float64,1}    # This needs to be a 2d array
    tstep::Float64
    time::DateTime
    
end


function NoSnow(tstep::Float64, time::DateTime, frac_lus::Array{Float64,1})
    
    p_in  = zeros(Float64, length(frac_lus))
    tair  = zeros(Float64, length(frac_lus))
    q_out = zeros(Float64, length(frac_lus))
    
    NoSnow(p_in, tair, q_out, frac_lus, tstep, time)
    
end


function get_param_ranges(model::NoSnow)
    
    param_range = Dict()
    
end


function init_states!(model::NoSnow, init_time::DateTime)
    
    model.time = init_time

    return nothing

end


function run_timestep(m::NoSnow)
    
    for reg in eachindex(m.frac_lus)
        
        m.q_out[reg] = m.p_in[reg]
        
    end
    
    m.time += Dates.Hour(m.tstep)
    
    return nothing
    
end
