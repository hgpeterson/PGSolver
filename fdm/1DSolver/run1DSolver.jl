using SparseArrays, LinearAlgebra, Printf, HDF5, PyPlot, PyCall

plt.style.use("~/paper_plots.mplstyle")
close("all")
pygui(false)

include("../../myJuliaLib.jl")
include("setParams1D.jl")
include("plottingLib1D.jl")
include("inversion1D.jl")
include("evolution1D.jl")
include("canonical1D.jl")

################################################################################
# run evolution integrations
################################################################################
#= b = evolve(500) =#
#= profilePlot(["b1000.h5", "b2000.h5", "b3000.h5", "b4000.h5", "b5000.h5"], 1) =#

#= sol = evolveCanonical1D(40000) =#
#= sol = steadyState() =#
profilePlot(["sol10.h5", "sol20.h5", "sol30.h5", "sol40.h5", "sol50.h5", "solSteady.h5"], 1)
