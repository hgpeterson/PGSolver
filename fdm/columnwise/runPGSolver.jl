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
# manually set `b` if you want
################################################################################
#= b = @. h*N^2*exp(-(z + H(x))/h) =#

#= file = h5open("/home/hpeter/Documents/ResearchCallies/rapid_adjustment/sims/sim000/b5000.h5") =#
#= b = read(file, "b") =#
#= close(file) =#

################################################################################
# test 1DAdjusted inversion
################################################################################
#= inversionLHS = getInversionLHS() =#
#= chi, uξ, uη, uσ, U = invert(b, inversionLHS) =#
#= u, v, w = transformFromTF(uξ, uη, uσ) =#
#= ridgePlot(chi, b, "chi", "") =#
#= savefig("chi.png") =#
#= close() =#
#= ridgePlot(u, b, "u", "") =#
#= savefig("u.png") =#
#= close() =#
#= ridgePlot(v, b, "v", "") =#
#= savefig("v.png") =#
#= close() =#
#= ridgePlot(w, b, "w", "") =#
#= savefig("w.png") =#
#= close() =#

#= inversionLHS1DAdjusted = getInversionLHS1DAdjusted() =#
#= inversionRHS1DAdjusted = getInversionRHS1DAdjusted(b) =#
#= sol = inversionLHS1DAdjusted\inversionRHS1DAdjusted =#
#= chi1D, u1D, v1D, w1D = postProcess1DAdjusted(sol) =#
#= ridgePlot(chi1D, b, "chi 1D", ""; vext=maximum(abs.(chi))) =#
#= savefig("chi1D_U.png") =#
#= close() =#
#= ridgePlot(u1D, b, "u 1D", ""; vext=maximum(abs.(u))) =#
#= savefig("u1D_U.png") =#
#= close() =#
#= ridgePlot(v1D, b, "v 1D", ""; vext=maximum(abs.(v))) =#
#= savefig("v1D_U.png") =#
#= close() =#
#= ridgePlot(w1D, b, "w 1D", ""; vext=maximum(abs.(w))) =#
#= savefig("w1D_U.png") =#
#= close() =#

################################################################################
# run evolution integrations
################################################################################

println("Computing inversion matrices")
inversionLHSs = Array{Any}(undef, nξ)
for i=1:nξ
    inversionLHSs[i] = lu(getInversionLHS(κ[i, :], H(ξ[i])))
end 

b = evolve(5000)

#= path = "" =#
#= dfiles = string.(path, ["b1000.h5", "b2000.h5", "b3000.h5", "b4000.h5", "b5000.h5"]) =#
#= profilePlot(dfiles, 1) =#


# tests
#= iξ = 150 =#
#= inversionLHS = getInversionLHS(κ[iξ, :], H(ξ[iξ])) =#
#= inversionRHS = getInversionRHS(b) =#


#= sol = zeros(nξ, nσ+1) =#
#= inversionRHS = getInversionRHS(b) =#
#= for i=1:nξ =#
#=     sol[i, :] = inversionLHSs[i]\inversionRHS[i, :] =#
#= end =#

#= sol = inversionLHS\inversionRHS[iξ, :] =#
#= println(size(inversionLHS)) =#
#= println(size(inversionRHS)) =#
#= sol = (inversionLHS\inversionRHS')' =#

#= sol = sol[iξ, :] =#
#= chi = sol[1:nσ] =#
#= U = sol[nσ+1] =#
#= plot(chi, σ) =#
#= axvline(U, lw=1, c="tab:orange") =#

#= chi, uξ, uη, uσ, U = postProcess(sol) =#
#= u, v, w = transformFromTF(uξ, uη, uσ) =#
#= ridgePlot(chi, b, "chi", "chi") =#
#= savefig("chi.png") =#
#= close() =#
#= ridgePlot(u, b, "u", "u") =#
#= savefig("u.png") =#
#= close() =#
#= ridgePlot(v, b, "v", "v") =#
#= savefig("v.png") =#
#= close() =#
#= ridgePlot(w, b, "w", "w") =#
#= savefig("w.png") =#
#= close() =#
