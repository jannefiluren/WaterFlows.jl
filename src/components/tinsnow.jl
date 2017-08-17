# Snow component

mutable struct TinSnow <: AbstractSnow
    
    swe::Array{Float64,2}
    tth::Float64
    ddf::Float64
    pcorr::Float64
    p_in::Array{Float64,1}
    tair::Array{Float64,1}
    q_out::Array{Float64,2}
    frac::Array{Float64,2}
    tstep::Float64
    time::DateTime
    
end


function TinSnow(tstep::Float64, time::DateTime, metadata::DataFrame)

    frac = convert(Array{Float64,2}, metadata)
    frac = transpose(frac)
    
    swe   = zeros(Float64, size(frac))
    p_in  = zeros(Float64, size(frac,2))
    tair  = zeros(Float64, size(frac,2))
    q_out = zeros(Float64, size(frac))
    
    tth = 0.0
    ddf = 3.69
    pcorr = 1.02
    
    TinSnow(swe, tth, ddf, pcorr, p_in, tair, q_out, frac, tstep, time)
    
end


function get_param_ranges(model::TinSnow)
    
    param_range = Dict(:tth => (-3.0, 3.0),
    :ddf => (0.1, 10.0),
    :pcorr => (0.5, 2.0))
    
end


function init_states!(model::TinSnow)

    for i in eachindex(model.swe)
        model.swe[i] = 0.0
    end
    
end


function run_timestep(m::TinSnow)
    
    for ireg in 1:size(m.frac, 2)

        # Compute solid and liquid precipitation
            
        psolid, pliquid = split_prec(m.p_in[ireg], m.tair[ireg], m.tth)
        
        psolid  = psolid * m.pcorr
        
        # Compute potential melt
            
        melt_pot = pot_melt(m.tair[ireg], m.ddf, m.tth)

        for ilus in 1:size(m.frac, 1)

            # Limit melt to available snow
            
            melt = min(m.swe[ilus, ireg], melt_pot)
            
            # Update snow water equivalents
            
            m.swe[ilus, ireg] += psolid - melt
            
            # Compute snowpack runoff
            
            m.q_out[ilus, ireg] = melt + pliquid
            
        end

    end
    
    m.time += Dates.Hour(m.tstep)
    
    return nothing
    
end



