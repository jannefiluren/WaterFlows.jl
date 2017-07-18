

mutable struct InputPTE <: AbstractInput

    prec::Array{Float64,2}
    tair::Array{Float64,2}
    epot::Array{Float64,1}

end


mutable struct InputPT <: AbstractInput

    prec::Array{Float64,2}
    tair::Array{Float64,2}

end

