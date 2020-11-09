using PyPlot, PyCall, Printf, SparseArrays, LinearAlgebra, HDF5

plt.style.use("~/paper_plots.mplstyle")
close("all")
pygui(false)

include("../../myJuliaLib.jl")
include("setParams.jl")
include("terrainFollowing.jl")
include("plottingLib.jl")
include("inversion.jl")
include("evolution.jl")

################################################################################
# run evolution integrations
################################################################################

b = evolve(5000; ξVariation=true)
#= profilePlot(["b1000.h5", "b2000.h5", "b3000.h5", "b4000.h5", "b5000.h5"], 1) =#
