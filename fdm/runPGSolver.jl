using PyPlot, PyCall, Printf, SparseArrays, LinearAlgebra, HDF5

plt.style.use("~/paper_plots.mplstyle")
close("all")
pygui(false)

include("../myJuliaLib.jl")
include("setParams.jl")
include("terrainFollowing.jl")
include("plottingLib.jl")
include("inversion.jl")
include("evolution.jl")

################################################################################
# run evolution integrations
################################################################################

#= print("Computing inversion matrices: ") =#
#= inversionLHSs = Array{Any}(undef, nξ) =#
#= for i=1:nξ =#
#=     inversionLHSs[i] = lu(getInversionLHS(κ[i, :], H(ξ[i]))) =#
#= end =# 
#= println("Done.") =#

#= b = evolve(5000) =#

################################################################################
# plots
################################################################################

#= path = "" =#
#= dfiles = string.(path, ["checkpoint1000.h5", "checkpoint2000.h5", "checkpoint3000.h5", "checkpoint4000.h5", "checkpoint5000.h5"]) =#
#= profilePlot(dfiles, 1) =#

include("talkPlots.jl")
#= folder = "/home/hpeter/Documents/ResearchCallies/rapid_adjustment/sims/sim008/" # bi κ =#
#= folder = "/home/hpeter/Documents/ResearchCallies/rapid_adjustment/sims/sim009/" # const κ =#
folder = "/home/hpeter/Documents/ResearchCallies/rapid_adjustment/sims/sim011/" # Pr
#= uAnimation(folder) =#
#= vAnimation(folder) =#
#= idealRidge() =#
#= uBalance(folder) =#
#= uvRidge(folder) =#
#= uvAnimation(folder) =#
#= uvPrScaling(folder) =#
#= BzChiPrScaling(folder) =#
pressureRidgePlots(string(folder, "full2D/Pr1/checkpoint1000.h5"))
