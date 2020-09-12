using PyPlot, HDF5

include("setParams.jl")
include("inversion.jl")
include("plottingLib.jl")

close("all")
plt.style.use("~/paper_plots.mplstyle")

# left-hand side for inversion equations
#= inversionLHS = lu(getInversionLHS()) =#

function compareChiEkman()
    file = h5open("b1000.h5", "r")
    b = read(file, "b")
    t = read(file, "t")
    close(file)

    # full buoyancy for isopycnals
    B = N^2*z + b 

    # chi and chiEkman
    chi, uξ, uη, uσ, U = invert(b, inversionLHS)
    chiEkman = getChiEkman(b)

    # plot both
    fig, ax = subplots(1, 3, figsize=(6.5, 6.5/1.62/1.8), gridspec_kw=Dict("width_ratios" =>[3, 3, 2]))

    vmax = maximum(abs.(chi))
    vmin = -vmax

    img = ax[1].pcolormesh(x/1000, z, chi, cmap="RdBu_r", vmin=vmin, vmax=vmax, rasterized=true)
    img = ax[2].pcolormesh(x/1000, z, chiEkman, cmap="RdBu_r", vmin=vmin, vmax=vmax, rasterized=true)
    fig.colorbar(img, ax=ax[1:2], label=L"$\chi$ (m$^2$ s$^{-1}$)", location="bottom", shrink=0.5, pad=0.2)
    #= fig.colorbar(img, ax=ax[2], label=L"$\chi$ (m$^2$ s$^{-1}$)") =#
    #= cax = fig.add_axes([0.01, 0.01, 6/7, 0.1]) =#
    #= fig.colorbar(img, cax=cax, label=L"$\chi$ (m$^2$ s$^{-1}$)") =#

    # isopycnal contours
    nLevels = 20
    lowerLevel = N^2*minimum(z)
    upperLevel = 0
    levels = lowerLevel:(upperLevel - lowerLevel)/(nLevels - 1):upperLevel
    ax[1].contour(x/1000, z, B, levels=levels, colors="k", alpha=0.3, linestyles="-")
    ax[2].contour(x/1000, z, B, levels=levels, colors="k", alpha=0.3, linestyles="-")

    # ridge shading
    ax[1].fill_between(x[:, 1]/1000, z[:, 1], minimum(z), color="k", alpha=0.3)
    ax[2].fill_between(x[:, 1]/1000, z[:, 1], minimum(z), color="k", alpha=0.3)

    # labels
    ax[1].set_title("numerical streamfunction")
    ax[1].set_xlabel(L"$x$ (km)")
    ax[1].set_ylabel(L"$z$ (m)")
    ax[2].set_title("analytical streamfunction")
    ax[2].set_xlabel(L"$x$ (km)")
    ax[3].set_xlabel(L"$\chi$ (m$^2$ s$^{-1}$)")

    # 1D plot
    ax[3].plot(chi[1, :], z[1, :], label="numerical")
    ax[3].plot(chiEkman[1, :], z[1, :], label="analytical", ls="--")
    ax[3].legend()

    #= tight_layout() =#

    savefig("chi_vs_chiEkman.pdf", bbox_inches="tight")
    #= savefig("chi_vs_chiEkman.pdf") =#
end

function compareRC20ridgePlots()
    # read
    file = h5open("b1000.h5", "r")
    b = read(file, "b")
    t = read(file, "t")
    close(file)

    # invert
    chi, uξ, uη, uσ, U = invert(b, inversionLHS)

    # transform
    u, v, w = transformFromTF(uξ, uη, uσ)

    # plot
    ridgePlot(u, b, "cross-ridge velocity", L"$u$ (m s$^{-1}$)"; vext=5e-5)
    savefig("u1000.pdf")
    close()
    ridgePlot(v, b, "along-ridge velocity", L"$v$ (m s$^{-1}$)"; vext=2e-2)
    savefig("v1000.pdf")
    close()
end

function compareRC20profilePlots()
    # read
    file = h5open("b1000.h5", "r")
    b1000 = read(file, "b")
    close(file)
    file = h5open("b2000.h5", "r")
    b2000 = read(file, "b")
    close(file)
    file = h5open("b3000.h5", "r")
    b3000 = read(file, "b")
    close(file)
    file = h5open("b4000.h5", "r")
    b4000 = read(file, "b")
    close(file)
    file = h5open("b5000.h5", "r")
    b5000 = read(file, "b")
    close(file)

    # invert
    chi1000, uξ1000, uη1000, uσ1000, U1000 = invert(b1000, inversionLHS)
    chi2000, uξ2000, uη2000, uσ2000, U2000 = invert(b2000, inversionLHS)
    chi3000, uξ3000, uη3000, uσ3000, U3000 = invert(b3000, inversionLHS)
    chi4000, uξ4000, uη4000, uσ4000, U4000 = invert(b4000, inversionLHS)
    chi5000, uξ5000, uη5000, uσ5000, U5000 = invert(b5000, inversionLHS)

    # transform
    u1000, v1000, w1000 = transformFromTF(uξ1000, uη1000, uσ1000)
    u2000, v2000, w2000 = transformFromTF(uξ2000, uη2000, uσ2000)
    u3000, v3000, w3000 = transformFromTF(uξ3000, uη3000, uσ3000)
    u4000, v4000, w4000 = transformFromTF(uξ4000, uη4000, uσ4000)
    u5000, v5000, w5000 = transformFromTF(uξ5000, uη5000, uσ5000)

    # strat
    bz1000 = zDerivativeTF(b1000)
    bz2000 = zDerivativeTF(b2000)
    bz3000 = zDerivativeTF(b3000)
    bz4000 = zDerivativeTF(b4000)
    bz5000 = zDerivativeTF(b5000)
    
    # index
    iξ = 1

    # plot strat
    fig, ax = plt.subplots(1)
    ax.plot(N^2 .+ bz1000[iξ, :], z[iξ, :], label="Day 1000")
    ax.plot(N^2 .+ bz2000[iξ, :], z[iξ, :], label="Day 2000")
    ax.plot(N^2 .+ bz3000[iξ, :], z[iξ, :], label="Day 3000")
    ax.plot(N^2 .+ bz4000[iξ, :], z[iξ, :], label="Day 4000")
    ax.plot(N^2 .+ bz5000[iξ, :], z[iξ, :], label="Day 5000")
    ax.legend(fontsize=6)
    ax.set_xlim([0, 1.4e-6])
    ax.set_xlabel(L"$B_z$, (s$^{-2}$)")
    ax.set_ylabel(L"$z$, (m)")
    ax.set_title("stratification")
    ax.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
    tight_layout()
    savefig("strat.pdf")
end

#= function profilePlotInit() =#
#=     fig, ax = subplots(2, 2, sharey=true) =#

#=     ax[1, 1].set_xlabel(L"$u$, (m s$^{-1}$)") =#
#=     ax[1, 1].set_ylabel(L"$z$, (m)") =#
#=     ax[1, 1].set_title("cross-ridge velocity") =#

#=     ax[1, 2].set_xlabel(L"$v$, (m s$^{-1}$)") =#
#=     #1= ax[1, 2].set_ylabel(L"$z$, (m)") =1# =#
#=     ax[1, 2].set_title("along-ridge velocity") =#

#=     ax[2, 1].set_xlabel(L"$w$, (m s$^{-1}$)") =#
#=     ax[2, 1].set_ylabel(L"$z$, (m)") =#
#=     ax[2, 1].set_title("vertical velocity") =#

#=     ax[2, 2].set_xlabel(L"$B_z$, (s$^{-2}$)") =#
#=     #1= ax[2, 2].set_ylabel(L"$z$, (m)") =1# =#
#=     ax[2, 2].set_title("stratification") =#

#=     tight_layout() =#

#=     ax[1, 1].ticklabel_format(style="sci", axis="x", scilimits=(0, 0)) =#
#=     ax[1, 2].ticklabel_format(style="sci", axis="x", scilimits=(0, 0)) =#
#=     ax[2, 1].ticklabel_format(style="sci", axis="x", scilimits=(0, 0)) =#
#=     ax[2, 2].ticklabel_format(style="sci", axis="x", scilimits=(0, 0)) =#
    
#=     return ax =#
#= end =#

#= function profilePlot(ax, uξ, uη, uσ, b, iξ, day) =#
#=     # convert to physical coordinates =# 
#=     u, v, w = transformFromTF(uξ, uη, uσ) =#

#=     # for stratification =#
#=     bz = zDerivativeTF(b) =#

#=     # plot =#
#=     ax[1, 1].plot(u[iξ, :], z[iξ, :], label=string("Day ", Int64(round(day)))) =#
#=     ax[1, 2].plot(v[iξ, :], z[iξ, :]) =#
#=     ax[2, 1].plot(w[iξ, :], z[iξ, :]) =#
#=     ax[2, 2].plot(N^2 .+ bz[iξ, :], z[iξ, :]) =#
#= end =#

#= """ =#
#=     plotCurrentState(t, chi, chiEkman, uξ, uη, uσ, b, iImg) =#

#= Plot the buoyancy and velocity state of the model at time `t` using label number `iImg`. =#
#= """ =#
#= function plotCurrentState(t, chi, chiEkman, uξ, uη, uσ, b, iImg) =#
#=     # convert to physical coordinates =# 
#=     u, v, w = transformFromTF(uξ, uη, uσ) =#

#=     # plots =#
#=     ridgePlot(chi, b, @sprintf("streamfunction at t = %.1f days", t/86400), L"$\chi$ (m$^2$ s$^{-1}$)") =#
#=     savefig(@sprintf("chi%03d.png", iImg)) =#
#=     close() =#

#=     ridgePlot(chiEkman, b, @sprintf("streamfunction theory at t = %.1f days", t/86400), L"$\chi$ (m$^2$ s$^{-1}$)") =#
#=     savefig(@sprintf("chiEkman%03d.png", iImg)) =#
#=     close() =#

#=     ridgePlot(b, b, @sprintf("buoyancy perturbation at t = %.1f days", t/86400), L"$b$ (m s$^{-2}$)") =#
#=     savefig(@sprintf("b%03d.png", iImg)) =#
#=     close() =#

#=     #1= ridgePlot(u, b, @sprintf("cross-ridge velocity at t = %.1f days", t/86400), L"$u$ (m s$^{-1}$)"; vext=5e-5) =1# =#
#=     ridgePlot(u, b, @sprintf("cross-ridge velocity at t = %.1f days", t/86400), L"$u$ (m s$^{-1}$)") =#
#=     savefig(@sprintf("u%03d.png", iImg)) =#
#=     close() =#

#=     #1= ridgePlot(v, b, @sprintf("along-ridge velocity at t = %.1f days", t/86400), L"$v$ (m s$^{-1}$)"; vext=2e-2) =1# =#
#=     ridgePlot(v, b, @sprintf("along-ridge velocity at t = %.1f days", t/86400), L"$v$ (m s$^{-1}$)") =#
#=     savefig(@sprintf("v%03d.png", iImg)) =#
#=     close() =#

#=     ridgePlot(w, b, @sprintf("vertical velocity at t = %.1f days", t/86400), L"$w$ (m s$^{-1}$)") =#
#=     savefig(@sprintf("w%03d.png", iImg)) =#
#=     close() =#
#= end =#


#= compareChiEkman() =#
#= compareRC20ridgePlots() =#
compareRC20profilePlots()
