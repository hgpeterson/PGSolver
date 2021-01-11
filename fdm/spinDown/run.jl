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

#= sol = evolve(5*τ) =#

################################################################################
# plots
################################################################################

#= path = "" =#
#= iSaves = 0:5 =#
#= dfiles = string.(path, "checkpoint", iSaves, ".h5") =#
#= profilePlot(dfiles) =#

#= MR91plot(dfiles) =#

#= folder = "/home/hpeter/ResearchCallies/rapid_adjustment/sims/sim014/" # small S =#
folder = "/home/hpeter/ResearchCallies/rapid_adjustment/sims/sim015/" # larger S
canonicalVsNoTransport(folder)
