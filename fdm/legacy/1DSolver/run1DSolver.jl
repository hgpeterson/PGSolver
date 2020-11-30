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

#= sol = evolveCanonical1D(40000) =#
#= sol = evolveCanonical1D(500) =#
sol = steadyState()
profilePlot(["sol1000.h5", "sol2000.h5", "sol3000.h5", "sol4000.h5", "sol5000.h5", "solSteady.h5"], 1)
