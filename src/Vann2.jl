module Vann2

using BlackBoxOptim
using CSV
using Distributions
using DataFrames

# Abstract types

abstract type AbstractModel end
abstract type AbstractSnow end
abstract type AbstractGlacier end
abstract type AbstractHydro end
abstract type AbstractInput end

export SemiDistComp, SemiDistFull
export SnowDistModel
export InputPTE, InputPT
export Gr4j, Hbv
export HbvLight
export TinSnow, NoSnow
export TinGlacier, HockGlacier, NoGlacier

export run_timestep, run_model
export load_data, crop_data
export get_param_ranges
export epot_zero, hamon, oudin
export init_states!
export run_model_calib
export get_params

include("inputs.jl")
include("components/gr4j.jl")
include("components/hbv.jl")
include("components/hbvlight.jl")
include("components/nosnow.jl")
include("components/tinsnow.jl")
include("components/noglacier.jl")
include("components/tinglacier.jl")
include("components/hockglacier.jl")
include("models/semidist.jl")
include("models/snowdist.jl")
include("utils_snow.jl")
include("solar_rad.jl")
include("utils_data.jl")
include("pot_evap.jl")
include("utils_calib.jl")

end
