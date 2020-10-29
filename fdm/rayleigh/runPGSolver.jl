include("setParams.jl")
include("plottingLib.jl")
include("inversion.jl")
include("evolution.jl")

################################################################################
# run evolution integrations
################################################################################

#= b = evolve(5000) =#
profilePlot(["b1000.h5", "b2000.h5", "b3000.h5", "b4000.h5", "b5000.h5"], 1)
