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

#= include("setParams.jl") =#
#= û, v, b, Px = evolve(5*τ_A) =#

#= τ_As = 10 .^range(1, 4, length=2^5) =#
#= τ_Ss = 10 .^range(1, 4, length=2^5) =#
#= Eks = 1 ./τ_Ss.^2 =#
#= Ss = 1 ./τ_As =#
#= vs = zeros((size(τ_Ss, 1), size(τ_As, 1))) =#

#= canonical = false =#
#= Pr = 1e3 =#
#= κ0 = 0 =#
#= κ1 = 1e-3 =#
#= h = 0 =#
#= v0 = -1 =#
#= bottomIntense = false =#
#= adaptiveTimestep = false =#
#= α = 0.5 =#
#= for i=1:size(τ_Ss, 1) =#
#=     println("i = ", i) =#
#=     for j=1:size(τ_As, 1) =#
#=         global Ek = Eks[i] =#
#=         global S = Ss[j] =#
#=         global τ_S = τ_Ss[i] =#
#=         global τ_A = τ_As[j] =#
#=         global Δt = minimum([τ_S/100, τ_A/100, 10]) =#
#=         global H = 1/sqrt(Ek) =#
#=         if H <= 100 =#
#=             global nẑ = 2^5 =#
#=         else =#
#=             global nẑ = 2^10 =#
#=         end =#
#=         global κ = κ1*ones(nẑ) =#
#=         global ẑ = @. H*(1 - cos(pi*(0:nẑ-1)/(nẑ-1)))/2 =#

#=         global tSave = 10*τ_A =#
#=         #1= global tSave = τ_A =1# =#
        
#=         û, v, b, Px = evolve(5*τ_A) =#
#=         vs[i, j] = v[end] =#
#=         #1= profilePlot(string.("checkpoint", 0:1:5, ".h5"); fname=string("profiles", i, j, ".png")) =1# =#
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
#= #1= img = ax.pcolormesh(τ_Ss, τ_As, vs'/v0, rasterized=true, shading="auto", vmin=0, vmax=1) =1# =#
#= #1= cb = colorbar(img, ax=ax, label=L"far-field $v/v_0$ at $t = 5\tau_A$") =1# =#
#= img = ax.contour(τ_Ss, τ_As, vs'/v0, levels=10) =#
#= cb = colorbar(img, ax=ax, label=L"far-field $v/v_0$ at $t = 5\tau_A$") =#
#= ax.loglog([0, 1], [0, 1], transform=ax.transAxes, "k--", lw=0.5) =#
#= tight_layout() =#
#= savefig("farfieldv.png") =#

################################################################################
# plots
################################################################################

#= path = "" =#
#= iSaves = 0:1:5 =#
#= dfiles = string.(path, "checkpoint", iSaves, ".h5") =#
#= profilePlot(dfiles) =#
