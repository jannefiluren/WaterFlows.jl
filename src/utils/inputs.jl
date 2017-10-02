

mutable struct InputPTE <: AbstractInput

    time::Array{DateTime,1}
    prec::Array{Float64,2}
    tair::Array{Float64,2}
    epot::Array{Float64,1}

end


mutable struct InputPT <: AbstractInput

    time::Array{DateTime,1}
    prec::Array{Float64,2}
    tair::Array{Float64,2}

end

