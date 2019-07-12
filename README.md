# WaterFlows

| **Build Status**                                                                                |
|:----------------------------------------------------------------------------------------------- |
 [![][travis-img]][travis-url] [![][codecov-img]][codecov-url]

Hydrological models implemented in the Julia programming language. The package currently provides a set of conceptual hydrological models including methods for computing potential evapotranspiration. The models are built using different components and can be combined in different constellations.

## Installation

The package can be installed with the Julia package manager. From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

````julia
pkg> add WaterFlows
````

And load the package using the command:

````julia
using WaterFlows
````

## Load input data

WaterFlows currently reads data in a specific text format (see examples for [Atnasjø](https://github.com/jmgnve/WaterFlows/tree/master/data/atnasjo) and [Fetvatn](https://github.com/jmgnve/WaterFlows/tree/master/data/fetvatn)). 

First, read the data for one of the example datasets, here Atnasjø:

````julia
path = joinpath(dirname(pathof(WaterFlows)), "..", "data", "atnasjo")
date, tair, prec, q_obs, frac_lus, frac_area, elev = load_data(path)
````

Second, compute potential evapotranspiration for the catchment:

````julia
lat = 60.0
epot = oudin(date, tair, lat, frac_area)
````

Third, create an object containing the input data that is required for running the models:

````julia
input = InputPTE(date, prec, tair, epot)
````

## Setting up and running a complete and predefined model

The most user-friendly method to run a model is to use a predefined model structure.

Start by specifying the time step length in hours and the time for the first input data: 

````julia
tstep = 24.0
tstart = date[1]
````

Next setup a model object containing necessary data required for running the model:

````julia
model = model_hbv_light(tstep, tstart, frac_lus)
````

Finally run the model:

````julia
q_sim = run_model(model, input)
````

## Building and running a model from components

A complete model can be built from existing components.

Again, start by specifying the time step length in hours and the time for the first input data: 

````julia
tstep = 24.0
tstart = date[1]
````

Next specify a snow, glacier and subsurface component:

````julia
snow = HbvLightSnow(tstep, tstart, frac_lus)
glacier = NoGlacier()
subsurf = Gr4j(tstep, tstart)
````

and create a model object:

````julia
model = ModelComp(snow, glacier, subsurf)
````

Finally run the model:

````julia
q_sim = run_model(model, input)
````

## Model calibration

A model can be calibrated by running:

````julia
param_tuned = run_model_calib(model, input, q_obs, warmup = 1, verbose = :verbose)
````

The model can be ran using the best-fit parameters with the following command:

````julia
set_params!(model, param_tuned)
q_sim = run_model(model, input)
````

## Available components

Currently, a complete model is splitted into three components representing snow melt, glacier melt and subsurface hydrological processes, respectivily.

Available snow components:

- HBV light snow melt component (HbvLightSnow)
- Simple temperature index snow melt component (TinSnow)
- Dummy component for neglecting snow melt (NoSnow)

Available glacier components:

- Radiation based glacier melt component (HockGlacier)
- Temperature index based glacier melt component (TinGlacier)
- Dummy component for neglecting glacier melt (NoGlacier)

Available subsurface components:

- A simple HBV subsurface model component (HbvLightSubsurf)
- HBV light subsurface component (Hbv)
- GR4J subsurface components (Gr4j)

The components can be ordered in any combination. However, note that the input arguments may differ between the components. For looking at the implementation of the components, click [here](https://github.com/jmgnve/WaterFlows/tree/master/src/components). The components are combined together to a complete model using [this code](https://github.com/jmgnve/WaterFlows/blob/master/src/models/model_components.jl).

## Available models

Currently only the HBV light model setup is available as a complete model as described above.


[travis-img]: https://travis-ci.org/jannefiluren/WaterFlows.jl.svg?branch=master
[travis-url]: https://travis-ci.org/jannefiluren/WaterFlows.jl

[codecov-img]: https://codecov.io/gh/jannefiluren/WaterFlows.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/jannefiluren/WaterFlows.jl
