# Examples

These scripts use [GLMakie](https://docs.makie.org/) for plotting, which is a
dependency of this example environment rather than of the `WaterFlows` package
itself. Set the environment up once from the repository root:

```julia
using Pkg
Pkg.activate("examples")
Pkg.develop(path = ".")   # add the local WaterFlows package
Pkg.instantiate()         # install GLMakie, DataFrames, Dates
```

Then run an example, e.g.:

```julia
include("examples/compare_different_timesteps.jl")
```
