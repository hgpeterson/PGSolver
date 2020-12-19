using PyPlot, PyCall, Printf, SparseArrays, LinearAlgebra, HDF5, Dierckx

plt.style.use("~/paper_plots.mplstyle")
close("all")
pygui(false)

include("../../myJuliaLib.jl")
include("setParams.jl")
include("evolution.jl")
include("plottingLib.jl")
include("rotated.jl")

################################################################################
# run evolution integrations
################################################################################

sol = evolve(10)

################################################################################
# plots
################################################################################

path = ""
#= tDays = 0:1000:5000 =#
#= tDays = 0:100:500 =#
#= tDays = 0:10:50 =#
tDays = 0:2:10
dfiles = string.(path, "checkpoint", tDays, ".h5")
profilePlot(dfiles)
