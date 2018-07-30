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

    iglacier = findall(names(frac_lus) .== :glacier)
    iglacier = iglacier[1]

    frac_lus = convert(Array{Float64,2}, frac_lus)
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
    
    for reg in eachindex(g.frac_lus)
        
        if s.swe[g.iglacier, 1] > 0.0
            g.q_out[reg] = 0.0
        else
            g.q_out[reg] = pot_melt(g.tair[reg], g.ddf)
        end
                
    end
    
    g.time += Dates.Hour(g.tstep)
    
    return nothing
    
end
