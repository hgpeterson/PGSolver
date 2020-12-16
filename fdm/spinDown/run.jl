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

sol = evolve(5000)

################################################################################
# plots
################################################################################

path = ""
dfiles = string.(path, ["checkpoint1000.h5", "checkpoint2000.h5", "checkpoint3000.h5", "checkpoint4000.h5", "checkpoint5000.h5"])
profilePlot(dfiles)
