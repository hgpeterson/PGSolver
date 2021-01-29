using PyPlot, PyCall, Printf, SparseArrays, LinearAlgebra, HDF5, Dierckx

plt.style.use("~/paper_plots.mplstyle")
close("all")
pygui(false)

include("../../../myJuliaLib.jl")
#= include("setParams.jl") =#
include("evolution.jl")
include("plottingLib.jl")
include("rotated.jl")

################################################################################
# run evolution integrations
################################################################################

include("setParams.jl")
sol = evolve(5*τ_A)

#= τ_As = 10 .^range(1, 4, length=2^7) =#
#= τ_Ss = 10 .^range(1, 4, length=2^7) =#
#= Eks = 1 ./τ_Ss.^2 =#
#= Ss = 1 ./τ_As =#
#= vs = zeros((size(τ_Ss, 1), size(τ_As, 1))) =#

#= canonical = false =#
#= Pr = 1e3 =#
#= κ0 = 0 =#
#= κ1 = 1e-3 =#
#= h = 0 =#
#= v0 = -1 =#
#= nẑ = 2^7 =#
#= bottomIntense = false =#
#= κ = κ1*ones(nẑ) =#
#= adaptiveTimestep = false =#
#= α = 0.5 =#
#= umap = reshape(1:3*nẑ, 3, nẑ) =#    
#= for i=1:size(τ_Ss, 1) =#
#=     println("i = ", i) =#
#=     for j=1:size(τ_As, 1) =#
#=         global Ek = Eks[i] =#
#=         global S = Ss[j] =#
#=         global τ_S = τ_Ss[i] =#
#=         global τ_A = τ_As[j] =#
#=         global Δt = τ_A =#
#=         global tSave = 10*τ_A =#
#=         global H = 1/sqrt(Ek) =#
#=         global ẑ = @. H*(1 - cos(pi*(0:nẑ-1)/(nẑ-1)))/2 =#

#=         sol = evolve(5*τ_A) =#
#=         v = sol[umap[2, :]] =#
#=         vs[i, j] = v[end] =#
#=     end =#
#= end =#

#= file = h5open("vs.h5", "w") =#
#= write(file, "vs", vs) =#
#= write(file, "τ_Ss", τ_Ss) =#
#= write(file, "τ_As", τ_As) =#
#= close(file) =#

#= file = h5open("vs.h5", "r") =#
#= vs = read(file, "vs") =#
#= τ_Ss = read(file, "τ_Ss") =#
#= τ_As = read(file, "τ_As") =#
#= close(file) =#

#= fig, ax = subplots(1) =#
#= ax.set_xlabel(L"$\tau_S$") =#
#= ax.set_ylabel(L"$\tau_A$") =#
#= img = ax.pcolormesh(τ_Ss, τ_As, vs'/v0, rasterized=true, shading="auto", vmin=0, vmax=1) =#
#= cb = colorbar(img, ax=ax, label=L"far-field $v/v_0$ at $t = 5\tau_A$") =#
#= #1= img = ax.contour(τ_Ss, τ_As, vs'/v0, levels=10) =1# =#
#= ax.loglog([0, 1], [0, 1], transform=ax.transAxes, "k--", lw=0.5) =#
#= tight_layout() =#
#= savefig("farfieldv.png") =#

################################################################################
# plots
################################################################################

path = ""
iSaves = 0:1:5
dfiles = string.(path, "checkpoint", iSaves, ".h5")
profilePlot(dfiles)
