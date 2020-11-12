using PyPlot, PyCall, Printf, SparseArrays, LinearAlgebra, HDF5, SpecialFunctions

plt.style.use("~/paper_plots.mplstyle")
close("all")
pygui(false)

include("../../../myJuliaLib.jl")
include("setParams.jl")
include("terrainFollowing.jl")
include("plottingLib.jl")
include("inversion.jl")
include("evolution.jl")

################################################################################
# run evolution integrations
################################################################################

#= b = evolve(5000) =#
#= profilePlot(["b1000.h5", "b2000.h5", "b3000.h5", "b4000.h5", "b5000.h5"], 1) =#
profilePlot(["b1000.h5", "b2000.h5", "b3000.h5"], 1)

for t=[1000 2000 3000 Inf]
    b1D, u1D, v1D, w1D = pointwise1D(t*86400)
    ridgePlot(b1D, b1D, "", "")
    savefig(string("b1D", t, ".png"))
end
