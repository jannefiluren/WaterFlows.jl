module WaterFlows

# External packages

using BlackBoxOptim
using CSV
using Distributions
using DataFrames
using Dates

# Abstract types

abstract type AbstractModel end
abstract type AbstractSnow end
abstract type AbstractGlacier end
abstract type AbstractSubsurf end
abstract type AbstractSubsurfLumped <: AbstractSubsurf end
abstract type AbstractSubsurfDist <: AbstractSubsurf end
abstract type AbstractInput end

# Input types

export InputPTE, InputPT

# Snow components

export TinSnow, HbvLightSnow, NoSnow

# Glacier components

export TinGlacier, HockGlacier, NoGlacier

# Subsurface components

export Gr4j, Hbv, HbvLightSubsurf

# Model setups

export ModelComp
export model_hbv_light
export model_gr4j

# Utilities

export run_timestep, run_model
export load_data, crop_data
export get_param_ranges, get_water_stored
export epot_zero, hamon, oudin
export init_states!, set_params!
export run_model_calib
export get_params
export nse, kge

# Utilities

include("utils/inputs.jl")
include("utils/utils_snow.jl")
include("utils/solar_rad.jl")
include("utils/utils_data.jl")
include("utils/pot_evap.jl")
include("utils/utils_calib.jl")

# Snow components

include("components/snow/nosnow.jl")
include("components/snow/tinsnow.jl")
include("components/snow/hbv_light_snow.jl")

# Glacier components

include("components/glacier/noglacier.jl")
include("components/glacier/tinglacier.jl")
include("components/glacier/hockglacier.jl")

# Subsurface components

include("components/subsurf/gr4j.jl")
include("components/subsurf/hbv.jl")
include("components/subsurf/hbv_light_subsurf.jl")

# Model setups

include("models/model_components.jl")
include("models/model_predefined.jl")

end
