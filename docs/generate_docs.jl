# Generate package documentation

using Weave

dir = Pkg.dir("Vann2")

cd(joinpath(dir, "docs"))

weave("README.jmd", out_path = dir, doctype = "github")
