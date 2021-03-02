using PyPlot, PyCall, Printf, SparseArrays, LinearAlgebra, HDF5, Dierckx, SpecialFunctions

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

print("Computing inversion matrices: ")
inversionLHSs = Array{Any}(undef, nξ)
for i=1:nξ
    inversionLHSs[i] = lu(getInversionLHS(κ[i, :], H(ξ[i])))
end 
# particular solution 
inversionRHS = getInversionRHS(zeros(nξ, nσ), 1)
solᵖ = computeSol(inversionRHS)
println("Done.")

b = evolve(5000)

################################################################################
# plots
################################################################################

path = ""
dfiles = string.(path, ["checkpoint1000.h5", "checkpoint2000.h5", "checkpoint3000.h5", "checkpoint4000.h5", "checkpoint5000.h5"])
profilePlot(dfiles, 1)
#= profilePlot(dfiles, Int64(round(nξ/2))) =#

#= include("talkPlots.jl") =#
#= folder = "/home/hpeter/Documents/ResearchCallies/rapid_adjustment/sims/sim008/" # bi κ =#
#= folder = "/home/hpeter/Documents/ResearchCallies/rapid_adjustment/sims/sim009/" # const κ =#
#= folder = "/home/hpeter/Documents/ResearchCallies/rapid_adjustment/sims/sim011/" # Pr (bi κ) =#
#= uvAnimation(folder) =#
#= chivAnimation(folder) =#
#= idealRidge() =#
#= uBalance(folder) =#
#= chiBalance(folder) =#
#= chiForSketch(folder) =#
#= ridge(folder) =#
#= uvPrScaling(folder) =#
#= BzChiPrScaling(folder) =#
#= pressureRidgePlots(string(folder, "full2D/Pr1/checkpoint1000.h5")) =#
#= profiles2Dvs1D(folder) =#
#= BzChi2DvsFixed(folder) =#
#= uv2DvsFixed(folder) =#
