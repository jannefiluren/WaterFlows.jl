# Glacier component

mutable struct HockGlacier <: AbstractGlacier
    
    f_m::Float64
    r_ice::Float64
    tair::Array{Float64,1}
    q_out::Array{Float64,1}
    frac_lus::Array{Float64,1}
    elev::Array{Float64,1}
    lat::Float64
    iglacier::Int64
    tstep::Float64
    time::DateTime
    
end


function HockGlacier(tstep::Float64, time::DateTime, frac_lus::DataFrame, lat::Float64, elev::Array{Float64,1})
    
    iglacier = findall(names(frac_lus) .== :glacier)
    iglacier = iglacier[1]
    
    frac_lus = convert(Array{Float64,2}, frac_lus)
    frac_lus = frac_lus[:, iglacier]
    
    tair  = zeros(Float64, length(frac_lus))
    q_out = zeros(Float64, length(frac_lus))
    
    f_m = 1.0
    r_ice = 1.0
    
    HockGlacier(f_m, r_ice, tair, q_out, frac_lus, elev, lat, iglacier, tstep, time)
    
end


function get_param_ranges(model::HockGlacier)
    
    param_range = Dict(
        :f_m => (0.1, 10.0),
        :r_ice => (0.1, 2.0))
    
end


function init_states!(model::HockGlacier, init_time::DateTime)

    model.time = init_time
    
    return nothing
    
end


function run_timestep(g::HockGlacier, s::AbstractSnow)
    
    for reg in eachindex(g.frac_lus)
        
        if s.swe[g.iglacier, 1] > 0.0
            g.q_out[reg] = 0.0
        elseif g.tair[reg] > 0.0
            R_so = clear_sky_rad(g.time, g.lat, g.elev[reg])
            g.q_out[reg] = (g.f_m + g.r_ice*R_so)*g.tair[reg]
        else
            g.q_out[reg] = 0.0
        end
        
    end
    
    g.time += Dates.Hour(g.tstep)
    
    return nothing
    
end
