################################################################################
# Functions useful for plotting
################################################################################

pl = pyimport("matplotlib.pylab")
inset_locator = pyimport("mpl_toolkits.axes_grid1.inset_locator")

"""
    ax = ridgePlot(field, b, titleString, cbarLabel; vext)

Create 2D plot of `field` with isopycnals given by the buoyancy perturbation `b`.
Set the title to `titleString` and colorbar label to `cbarLabel`. Return the axis 
handle `ax`.

Optional: set the vmin/vmax manually with vext.
"""
function ridgePlot(field, b, titleString, cbarLabel; vext=nothing, cmap="RdBu_r")
    # full buoyancy for isopycnals
    B = N^2*z + b 

    fig, ax = subplots(1)

    # set min and max
    if vext == nothing
        vmax = maximum(abs.(field))
        vmin = -vmax
        extend = "neither"
    else
        vmax = vext
        vmin = -vext
        extend = "both"
    end

    # regular min and max for viridis
    if cmap == "viridis"
        vmin = minimum(field)
        vmax = maximum(field)
        extend = "neither"
    end

    # 2D plot
    img = ax.pcolormesh(x/1000, z, field, cmap=cmap, vmin=vmin, vmax=vmax, rasterized=true)
    cb = colorbar(img, ax=ax, label=cbarLabel, extend=extend)
    cb.ax.ticklabel_format(style="sci", scilimits=(-3, 3))

    # isopycnal contours
    nLevels = 20
    lowerLevel = N^2*minimum(z)
    upperLevel = 0
    levels = lowerLevel:(upperLevel - lowerLevel)/(nLevels - 1):upperLevel
    ax.contour(x/1000, z, B, levels=levels, colors="k", alpha=0.3, linestyles="-", linewidths=0.5)

    # ridge shading
    ax.fill_between(x[:, 1]/1000, z[:, 1], minimum(z), color="k", alpha=0.3, lw=0.0)

    # labels
    ax.set_title(titleString)
    ax.set_xlabel(L"$x$ (km)")
    ax.set_ylabel(L"$z$ (m)")

    tight_layout()
    
    return ax
end

#= """ =#
#=     ax = ridgePlotNoIsopycnals(field, titleString, cbarLabel; vext) =#

#= Create 2D plot of `field` without isopycnals. =#
#= Set the title to `titleString` and colorbar label to `cbarLabel`. Return the axis =# 
#= handle `ax`. =#

#= Optional: set the vmin/vmax manually with vext. =#
#= """ =#
#= function ridgePlotNoIsopycnals(field, titleString, cbarLabel; vext=nothing, cmap="RdBu_r") =#
#=     fig, ax = subplots(1) =#

#=     # set min and max =#
#=     if vext == nothing =#
#=         vmax = maximum(abs.(field)) =#
#=         vmin = -vmax =#
#=         extend = "neither" =#
#=     else =#
#=         vmax = vext =#
#=         vmin = -vext =#
#=         extend = "both" =#
#=     end =#

#=     # regular min and max for viridis =#
#=     if cmap == "viridis" =#
#=         vmin = minimum(field) =#
#=         vmax = maximum(field) =#
#=         extend = "neither" =#
#=     end =#

#=     # 2D plot =#
#=     img = ax.pcolormesh(x/1000, z, field, cmap=cmap, vmin=vmin, vmax=vmax, rasterized=true) =#
#=     cb = colorbar(img, ax=ax, label=cbarLabel, extend=extend) =#
#=     cb.ax.ticklabel_format(style="sci", scilimits=(-3, 3)) =#

#=     # ridge shading =#
#=     ax.fill_between(x[:, 1]/1000, z[:, 1], minimum(z), color="k", alpha=0.3) =#

#=     # labels =#
#=     ax.set_title(titleString) =#
#=     ax.set_xlabel(L"$x$ (km)") =#
#=     ax.set_ylabel(L"$z$ (m)") =#

#=     tight_layout() =#
    
#=     return ax =#
#= end =#

#= """ =#
#=     ax = ridgePlot2Fields(field1, field2, titleString, cbarLabel; vext) =#

#= Create 2D plot of `field1` and `field2`. =#
#= Set the title to `titleString` and colorbar label to `cbarLabel`. Return the axis =# 
#= handle `ax`. =#

#= Optional: set the vmin/vmax manually with vext. =#
#= """ =#
#= function ridgePlot2Fields(field1, field2, titleString, cbarLabel; vext=nothing, cmap="RdBu_r") =#
#=     fig, ax = subplots(1) =#

#=     # set min and max =#
#=     if vext == nothing =#
#=         vmax = maximum(abs.(field1)) =#
#=         vmin = -vmax =#
#=         extend = "neither" =#
#=     else =#
#=         vmax = vext =#
#=         vmin = -vext =#
#=         extend = "both" =#
#=     end =#

#=     # regular min and max for viridis =#
#=     if cmap == "viridis" =#
#=         vmin = minimum(field1) =#
#=         vmax = maximum(field1) =#
#=         extend = "neither" =#
#=     end =#

#=     # 2D plot =#
#=     img = ax.pcolormesh(x/1000, z, field1, cmap=cmap, vmin=vmin, vmax=vmax, rasterized=true) =#
#=     cb = colorbar(img, ax=ax, label=cbarLabel, extend=extend) =#
#=     cb.ax.ticklabel_format(style="sci", scilimits=(-3, 3)) =#

#=     # contours =#
#=     ax.contour(x/1000, z, field2, levels=10, colors="k", alpha=0.3, linestyles="-") =#

#=     # ridge shading =#
#=     ax.fill_between(x[:, 1]/1000, z[:, 1], minimum(z), color="k", alpha=0.3) =#

#=     # labels =#
#=     ax.set_title(titleString) =#
#=     ax.set_xlabel(L"$x$ (km)") =#
#=     ax.set_ylabel(L"$z$ (m)") =#

#=     tight_layout() =#
    
#=     return ax =#
#= end =#

"""
    profilePlot(datafiles, iξ)

Plot profiles of b, u, v, w from HDF5 snapshot files of buoyancy in the `datafiles` list
at ξ = ξ[iξ].
"""
function profilePlot(datafiles, iξ)
    # init plot
    fig, ax = subplots(2, 2, figsize=(6.5, 6.5/1.62))

    # insets
    axins21 = inset_locator.inset_axes(ax[2, 1], width="40%", height="40%")
    axins22 = inset_locator.inset_axes(ax[2, 2], width="40%", height="40%")

    ax[1, 1].set_xlabel(L"$v$ (m s$^{-1}$)")
    ax[1, 1].set_ylabel(L"$z$ (m)")
    ax[1, 1].set_title("along-ridge velocity")

    ax[1, 2].set_xlabel(L"$B_z$ (s$^{-2}$)")
    ax[1, 2].set_ylabel(L"$z$ (m)")
    ax[1, 2].set_title("stratification")

    ax[2, 1].set_xlabel(L"$u$ (m s$^{-1}$)")
    ax[2, 1].set_ylabel(L"$z$ (m)")
    ax[2, 1].set_title("cross-ridge velocity")

    ax[2, 2].set_xlabel(L"$w$ (m s$^{-1}$)")
    ax[2, 2].set_ylabel(L"$z$ (m)")
    ax[2, 2].set_title("vertical velocity")

    #= tight_layout() =#
    subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.9, wspace=0.3, hspace=0.6)

    for a in ax
        a.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
    end
    axins21.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
    axins22.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))

    # left-hand side for inversion equations
    inversionLHS = lu(getInversionLHS())

    # color map
    colors = pl.cm.viridis(range(1, 0, length=5))

    # zoomed z
    ax[2, 1].set_ylim([z[iξ, 1], z[iξ, 1] + H0/10])
    ax[2, 2].set_ylim([z[iξ, 1], z[iξ, 1] + H0/10])

    # plot data from `datafiles
    for i=1:size(datafiles, 1)
        file = h5open(datafiles[i], "r")
        b = read(file, "b")
        t = read(file, "t")
        close(file)

        # invert buoyancy for flow
        chi, uξ, uη, uσ, U = invert(b, inversionLHS)

        # convert to physical coordinates 
        u, v, w = transformFromTF(uξ, uη, uσ)

        # stratification
        Bz = N^2 .+ zDerivativeTF(b)

        # colors and labels
        label = string("Day ", Int64(round(t/86400)))
        c = colors[i, :]

        # plot
        ax[1, 1].plot(v[iξ, :],  z[iξ, :], c=c, label=label)
        ax[1, 2].plot(Bz[iξ, :], z[iξ, :], c=c)
        ax[2, 1].plot(u[iξ, :],  z[iξ, :], c=c)
        axins21.plot(u[iξ, :],  z[iξ, :], c=c)
        ax[2, 2].plot(w[iξ, :],  z[iξ, :], c=c)
        axins22.plot(w[iξ, :],  z[iξ, :], c=c)
    end

    ax[1, 1].legend()

    savefig("profiles.png", bbox="inches")
end

#= """ =#
#=     profilePlot1DAdjusted(datafiles, iξ) =#

#= Plot profiles of b, u, v, w from HDF5 snapshot files of buoyancy in the `datafiles` list =#
#= at ξ = ξ[iξ]. =#
#= """ =#
#= function profilePlot1DAdjusted(datafiles, iξ) =#
#=     # init plot =#
#=     fig, ax = subplots(2, 2, figsize=(6.5, 6.5/1.62), sharey=true) =#

#=     ax[1, 1].set_xlabel(L"$u$ (m s$^{-1}$)") =#
#=     ax[1, 1].set_ylabel(L"$z$ (m)") =#
#=     ax[1, 1].set_title("cross-ridge velocity") =#

#=     ax[1, 2].set_xlabel(L"$v$ (m s$^{-1}$)") =#
#=     ax[1, 2].set_title("along-ridge velocity") =#

#=     ax[2, 1].set_xlabel(L"$w$ (m s$^{-1}$)") =#
#=     ax[2, 1].set_ylabel(L"$z$ (m)") =#
#=     ax[2, 1].set_title("vertical velocity") =#

#=     ax[2, 2].set_xlabel(L"$B_z$ (s$^{-2}$)") =#
#=     ax[2, 2].set_title("stratification") =#

#=     tight_layout() =#

#=     ax[1, 1].ticklabel_format(style="sci", axis="x", scilimits=(0, 0)) =#
#=     ax[1, 2].ticklabel_format(style="sci", axis="x", scilimits=(0, 0)) =#
#=     ax[2, 1].ticklabel_format(style="sci", axis="x", scilimits=(0, 0)) =#
#=     ax[2, 2].ticklabel_format(style="sci", axis="x", scilimits=(0, 0)) =#

#=     # left-hand side for inversion equations =#
#=     inversionLHS = lu(getInversionLHS1DAdjusted()) =#

#=     # color map =#
#=     colors = pl.cm.viridis(range(1, 0, length=5)) =#

#=     # plot data from `datafiles` =#
#=     for bfile in datafiles =#
#=         file = h5open(bfile, "r") =#
#=         b = read(file, "b") =#
#=         t = read(file, "t") =#
#=         close(file) =#

#=         # invert buoyancy for flow =#
#=         chi, uξ, uη, uσ, U = invert1DAdjusted(b, inversionLHS) =#

#=         # convert to physical coordinates =# 
#=         u, v, w = transformFromTF(uξ, uη, uσ) =#

#=         # stratification =#
#=         Bz = N^2 .+ zDerivativeTF(b) =#

#=         # colors and labels =#
#=         label = string("Day ", Int64(round(t/86400))) =#
#=         c = colors[i, :] =#

#=         # plot =#
#=         ax[1, 1].plot(u[iξ, :],  z[iξ, :], c=c, label=label) =#
#=         ax[1, 2].plot(v[iξ, :],  z[iξ, :], c=c) =#
#=         ax[2, 1].plot(w[iξ, :],  z[iξ, :], c=c) =#
#=         ax[2, 2].plot(Bz[iξ, :], z[iξ, :], c=c) =#
#=     end =#

#=     ax[1, 1].legend() =#

#=     fig.savefig("profiles.png") =#
#= end =#

#= """ =#
#=     profilePlot1D(datafiles, iξ) =#

#= Plot profiles of b, u, v, w from HDF5 snapshot files of buoyancy in the `datafiles` list =#
#= at ξ = ξ[iξ]. =#
#= """ =#
#= function profilePlot1D(datafiles, iξ) =#
#=     # init plot =#
#=     #1= fig, ax = subplots(2, 2, figsize=(6.5, 6.5/1.62), sharey=true) =1# =#
#=     fig, ax = subplots(2, 2, figsize=(6.5, 6.5/1.62)) =#

#=     ax[1, 1].set_xlabel(L"$u$ (m s$^{-1}$)") =#
#=     ax[1, 1].set_ylabel(L"$z$ (m)") =#
#=     ax[1, 1].set_title("cross-ridge velocity") =#

#=     ax[1, 2].set_xlabel(L"$v$ (m s$^{-1}$)") =#
#=     ax[1, 2].set_title("along-ridge velocity") =#

#=     ax[2, 1].set_xlabel(L"$w$ (m s$^{-1}$)") =#
#=     ax[2, 1].set_ylabel(L"$z$ (m)") =#
#=     ax[2, 1].set_title("vertical velocity") =#

#=     ax[2, 2].set_xlabel(L"$B_z$ (s$^{-2}$)") =#
#=     ax[2, 2].set_title("stratification") =#

#=     tight_layout() =#

#=     ax[1, 1].ticklabel_format(style="sci", axis="x", scilimits=(0, 0)) =#
#=     ax[1, 2].ticklabel_format(style="sci", axis="x", scilimits=(0, 0)) =#
#=     ax[2, 1].ticklabel_format(style="sci", axis="x", scilimits=(0, 0)) =#
#=     ax[2, 2].ticklabel_format(style="sci", axis="x", scilimits=(0, 0)) =#

#=     # color map =#
#=     colors = pl.cm.viridis(range(1, 0, length=5)) =#

#=     # plot data from `datafiles` =#
#=     for i=1:size(datafiles, 1) =#
#=         file = h5open(datafiles[i], "r") =#
#=         u = read(file, "u") =#
#=         v = read(file, "v") =#
#=         w = read(file, "w") =#
#=         b = read(file, "b") =#
#=         t = read(file, "t") =#
#=         close(file) =#

#=         # for stratification =#
#=         #1= bz = zDerivativeTF(b) =1# =#
#=         #1= Bz = N^2 .+ bz =1# =#
#=         bz = differentiate(b, z[iξ, :]) =#
#=         Bz = N^2*cosθ[iξ, 1] .+ bz =#

#=         # plot =#
#=         if t == Inf =#
#=             label = "Steady State" =#
#=             c = "k" =#
#=         else =#
#=             label = string("Day ", Int64(round(t/86400))) =#
#=             c = colors[i, :] =#
#=         end =#
#=         #1= ax[1, 1].plot(u[iξ, :], z[iξ, :], label=string("Day ", Int64(round(t/86400)))) =1# =#
#=         #1= ax[1, 2].plot(v[iξ, :], z[iξ, :]) =1# =#
#=         #1= ax[2, 1].plot(w[iξ, :], z[iξ, :]) =1# =#
#=         #1= ax[2, 2].plot(N^2 .+ bz[iξ, :], z[iξ, :]) =1# =#
#=         ax[1, 1].plot(u,  z[iξ, :], c=c, label=label) =#
#=         ax[1, 2].plot(v,  z[iξ, :], c=c) =#
#=         ax[2, 1].plot(w,  z[iξ, :], c=c) =#
#=         ax[2, 2].plot(Bz, z[iξ, :], c=c) =#
#=     end =#

#=     ax[1, 1].legend() =#
#=     ax[1, 1].set_ylim([-H(ξ[iξ]), -H(ξ[iξ]) + 200]) =#
#=     ax[1, 2].set_ylim([-H(ξ[iξ]), -H(ξ[iξ]) + 1000]) =#
#=     ax[2, 1].set_ylim([-H(ξ[iξ]), -H(ξ[iξ]) + 200]) =#
#=     ax[2, 2].set_ylim([-H(ξ[iξ]), -H(ξ[iξ]) + 1000]) =#

#=     fig.savefig("profiles.png") =#
#= end =#

"""
    plotCurrentState(t, chi, chiEkman, uξ, uη, uσ, b, iImg)

Plot the buoyancy and velocity state of the model at time `t` using label number `iImg`.
"""
function plotCurrentState(t, chi, chiEkman, uξ, uη, uσ, b, iImg)
    # convert to physical coordinates 
    u, v, w = transformFromTF(uξ, uη, uσ)

    # plots
    ridgePlot(chi, b, @sprintf("streamfunction at t = %4d days", t/86400), L"$\chi$ (m$^2$ s$^{-1}$)")
    savefig(@sprintf("chi%03d.png", iImg))
    close()

    ridgePlot(chiEkman, b, @sprintf("streamfunction theory at t = %4d days", t/86400), L"$\chi$ (m$^2$ s$^{-1}$)")
    savefig(@sprintf("chiEkman%03d.png", iImg))
    close()

    ridgePlot(b, b, @sprintf("buoyancy perturbation at t = %4d days", t/86400), L"$b$ (m s$^{-2}$)")
    savefig(@sprintf("b%03d.png", iImg))
    close()

    ridgePlot(u, b, @sprintf("cross-ridge velocity at t = %4d days", t/86400), L"$u$ (m s$^{-1}$)")
    savefig(@sprintf("u%03d.png", iImg))
    close()

    ridgePlot(v, b, @sprintf("along-ridge velocity at t = %4d days", t/86400), L"$v$ (m s$^{-1}$)")
    savefig(@sprintf("v%03d.png", iImg))
    close()

    ridgePlot(w, b, @sprintf("vertical velocity at t = %4d days", t/86400), L"$w$ (m s$^{-1}$)")
    savefig(@sprintf("w%03d.png", iImg))
    close()
end
#= function plotCurrentState(t, u, v, w, b, iImg) =#
#=     # plots =#
#=     ridgePlot(b, b, @sprintf("buoyancy perturbation at t = %4d days", t/86400), L"$b$ (m s$^{-2}$)") =#
#=     savefig(@sprintf("b%03d.png", iImg)) =#
#=     close() =#

#=     ridgePlot(u, b, @sprintf("cross-ridge velocity at t = %4d days", t/86400), L"$u$ (m s$^{-1}$)") =#
#=     savefig(@sprintf("u%03d.png", iImg)) =#
#=     close() =#

#=     ridgePlot(v, b, @sprintf("along-ridge velocity at t = %4d days", t/86400), L"$v$ (m s$^{-1}$)") =#
#=     savefig(@sprintf("v%03d.png", iImg)) =#
#=     close() =#

#=     ridgePlot(w, b, @sprintf("vertical velocity at t = %4d days", t/86400), L"$w$ (m s$^{-1}$)") =#
#=     savefig(@sprintf("w%03d.png", iImg)) =#
#=     close() =#
#= end =#
