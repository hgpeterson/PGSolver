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
#= iSaves = 0:1:5 =#
#= dfiles = string.(path, "checkpoint", iSaves, ".h5") =#
#= profilePlot(dfiles) =#

#= MR91plot(dfiles) =#

#= folder = "/home/hpeter/Documents/ResearchCallies/sims/sim014/" # small S =#
#= folder = "/home/hpeter/Documents/ResearchCallies/sims/sim015/" # larger S =#
#= folder = "/home/hpeter/Documents/ResearchCallies/sims/sim017/" # small S double H =#
folder = "/home/hpeter/Documents/ResearchCallies/sims/sim018/" # larger S double H
canonicalVsNoTransport(folder)
