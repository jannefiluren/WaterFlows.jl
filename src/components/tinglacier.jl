# Glacier component

mutable struct TinGlacier <: AbstractGlacier

    ddf::Float64
    tair::Array{Float64,1}
    q_out::Array{Float64,1}
    frac::Array{Float64,1}
    iglacier::Int64
    tstep::Float64
    time::DateTime
    
end


function TinGlacier(tstep::Float64, time::DateTime, metadata::DataFrame)

    iglacier = find(names(metadata) .== :glacier)
    iglacier = iglacier[1]

    frac = convert(Array{Float64,2}, metadata)
    frac = frac[:, iglacier]
    
    tair  = zeros(Float64, length(frac))
    q_out = zeros(Float64, length(frac))

    ddf = 5.0
    
    TinGlacier(ddf, tair, q_out, frac, iglacier, tstep, time)
    
end


function get_param_ranges(model::TinGlacier)
    
    param_range = Dict(:ddf => (0.1, 10.0))
    
end


function init_states!(model::TinGlacier)

    return nothing

end


function run_timestep(g::TinGlacier, s::AbstractSnow)
    
    for reg in eachindex(g.frac)
        
        if s.swe[g.iglacier, 1] > 0.0
            g.q_out[reg] = 0.0
        else
            g.q_out[reg] = pot_melt(g.tair[reg], g.ddf)
        end
                
    end
    
    g.time += Dates.Hour(g.tstep)
    
    return nothing
    
end
