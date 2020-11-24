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
# manually set `b` if you want
################################################################################
#= b = zeros(nξ, nσ) =#
#= bx = zeros(nξ, nσ) =#
#= for i=1:nξ =#
#=     for j=1:nσ =#
#=         # exponential in σ (centered at bottom) =#
#=         decay_scale = 200/H0 =#
#=         A = decay_scale*N^2*H0 =#
#=         b[i, j] = A*exp(-N^2*H(ξ[i])*(σ[j] + 1)/A) =#
#=         bx[i, j] = -N^2*Hx(ξ[i])*exp(-N^2*H(ξ[i])*(σ[j] + 1)/A) =#

#=         #1= # gaussian in σ (centered at bottom) =1# =#
#=         #1= b[i, j] = N^2*amp*exp(-(σ[j] + 1)^2/2/(1/4)^2) =1# =#
#=         #1= bx[i, j] = N^2*amp*exp(-(σ[j] + 1)^2/2/(1/4)^2)*(σ[j] + 1)/(1/4)^2*σ[j]*Hx(ξ[i])/H(ξ[i]) =1# =#
#=     end =#
#= end =#

#= b = @. h*N^2*exp(-(z + H(x))/h) =#
#= #1= bx = @. -N^2*Hx(x)*exp(-(z + H(x))/h) =1# =#
#= bx = zeros(nξ, nσ) =#
#= for j=1:nσ =#
#=     bx[:, j] = differentiate(b[:, j], ξ) =#
#= end =#
#= for i=1:nξ =#
#=     bx[i, :] .-= Hx(ξ[i])*σ.*differentiate(b[i, :], σ)/H(ξ[i]) =#
#= end =#

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

#= b = evolve(5000) =#

#= path = "" =#
#= dfiles = string.(path, ["b1000.h5", "b2000.h5", "b3000.h5", "b4000.h5", "b5000.h5"]) =#
#= profilePlot(dfiles, 1) =#

################################################################################
# plots
################################################################################
# advection terms
inversionLHS = lu(getInversionLHS())
figP, axP = subplots(2, 3, figsize=(6.5, 6.5/1.8), sharey=true)
folder = "../../sims/sim006/"
vmax = 0
pl = pyimport("matplotlib.pylab")
colors = pl.cm.viridis(range(1, 0, length=5))
bfiles = string.(folder, ["b1000.h5", "b2000.h5", "b3000.h5", "b4000.h5"])
for i=1:size(bfiles, 1)
    file = h5open(bfiles[i], "r")
    b = read(file, "b")
    t = read(file, "t")
    close(file)

    chi, uξ, uη, uσ, U = invert(b, inversionLHS) 

    adv1 = -uξ.*ξDerivativeTF(b)
    adv2 = -uσ.*σDerivativeTF(b)
    adv3 = -N^2*uξ.*Hx.(ξξ).*σσ
    adv4 = -N^2*uσ.*H.(ξξ)
    diff = σDerivativeTF(κ.*(N^2 .+ σDerivativeTF(b)./H.(ξξ)))./H.(ξξ)
    sum = adv1 + adv2 + adv3 + adv4 + diff

    vmax = maximum([maximum(abs.(adv1)) maximum(abs.(adv2)) maximum(abs.(adv3)) maximum(abs.(adv4)) maximum(abs.(diff))])

    fig, ax = subplots(2, 3, figsize=(6.5, 6.5/1.8), sharey=true, sharex=true)
    img = ax[1, 1].pcolormesh(ξξ/1000, σσ, adv1, vmin=-vmax, vmax=vmax, cmap="RdBu_r")
    colorbar(img, ax=ax[1, 1], label=L"-u^\xi b_\xi")
    img = ax[1, 2].pcolormesh(ξξ/1000, σσ, adv2, vmin=-vmax, vmax=vmax, cmap="RdBu_r")
    colorbar(img, ax=ax[1, 2], label=L"-u^\sigma b_\sigma")
    img = ax[1, 3].pcolormesh(ξξ/1000, σσ, adv3, vmin=-vmax, vmax=vmax, cmap="RdBu_r")
    colorbar(img, ax=ax[1, 3], label=L"-N^2u^\xi H_x\sigma")
    img = ax[2, 1].pcolormesh(ξξ/1000, σσ, adv4, vmin=-vmax, vmax=vmax, cmap="RdBu_r")
    colorbar(img, ax=ax[2, 1], label=L"-N^2u^\sigma H")
    img = ax[2, 2].pcolormesh(ξξ/1000, σσ, diff, vmin=-vmax, vmax=vmax, cmap="RdBu_r")
    cb = colorbar(img, ax=ax[2, 2], label=L"H^{-1}[\kappa(N^2 + H^{-1}b_\sigma)]_\sigma")
    cb.ax.ticklabel_format(style="sci", scilimits=(-3, 3))
    img = ax[2, 3].pcolormesh(ξξ/1000, σσ, sum, vmin=-vmax, vmax=vmax, cmap="RdBu_r")
    cb = colorbar(img, ax=ax[2, 3], label=L"b_t")
    cb.ax.ticklabel_format(style="sci", scilimits=(-3, 3))
    ax[1, 1].set_ylabel(L"\sigma")
    ax[2, 1].set_ylabel(L"\sigma")
    ax[2, 1].set_xlabel(L"\xi")
    ax[2, 2].set_xlabel(L"\xi")
    ax[2, 3].set_xlabel(L"\xi")
    fig.tight_layout()
    fig.savefig(@sprintf("evol%d.png", t/86400))

    axP[1, 1].plot(adv1[1, :], σ, c=colors[i, :], label=string("Day ", Int64(t/86400)))
    axP[1, 1].set_xlabel(L"-u^\xi b_\xi")
    axP[1, 2].plot(adv2[1, :], σ, c=colors[i, :])
    axP[1, 2].set_xlabel(L"-u^\sigma b_\sigma")
    axP[1, 3].plot(adv3[1, :], σ, c=colors[i, :])
    axP[1, 3].set_xlabel(L"-N^2u^\xi H_x\sigma")
    axP[2, 1].plot(adv4[1, :], σ, c=colors[i, :])
    axP[2, 1].set_xlabel(L"-N^2u^\sigma H")
    axP[2, 2].plot(diff[1, :], σ, c=colors[i, :])
    axP[2, 2].set_xlabel(L"H^{-1}[\kappa(N^2 + H^{-1}b_\sigma)]_\sigma")
    axP[2, 3].plot(sum[1, :], σ, c=colors[i, :])
    axP[2, 3].plot(adv3[1, :] + adv4[1, :] + diff[1, :], σ, c="k", ls=":")
    axP[2, 3].plot(adv3[1, :] + diff[1, :], σ, c=colors[i, :], ls="--")
    axP[2, 3].set_xlabel(L"b_t")
end
axP[1, 1].set_ylabel(L"\sigma")
axP[2, 1].set_ylabel(L"\sigma")
axP[2, 2].ticklabel_format(style="sci", scilimits=(0, 0))
axP[2, 3].ticklabel_format(style="sci", scilimits=(0, 0))
axP[1, 1].legend()
figP.tight_layout()
figP.savefig("evolProfiles_zoom.png")
for ax in axP
    ax.set_xlim([-3e-12, 3e-12])
end
figP.savefig("evolProfiles.png")
