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

    # no spines
    ax.spines["left"].set_visible(false)
    ax.spines["bottom"].set_visible(false)

    tight_layout()
    
    return ax
end

"""
    profilePlot(datafiles, ix)

Plot profiles of b, u, v, w from HDF5 snapshot files of buoyancy in the `datafiles` list
at ξ = ξ[ix].
"""
function profilePlot(datafiles, ix)
    # init plot
    fig, ax = subplots(2, 3, figsize=(10, 6.5/1.62))

    # insets
    axins21 = inset_locator.inset_axes(ax[2, 1], width="40%", height="40%")

    ax[1, 1].set_xlabel(L"$B_z$ (s$^{-2}$)")
    ax[1, 1].set_ylabel(L"$z$ (m)")
    ax[1, 1].set_title("stratification")

    ax[1, 2].set_xlabel(L"$\chi$ (m$^2$ s$^{-1}$)")
    ax[1, 2].set_ylabel(L"$z$ (m)")
    ax[1, 2].set_title("streamfunction")

    ax[1, 3].set_xlabel(L"$-b_{\hat{z}}\sin\theta$ (s$^{-2}$)")
    ax[1, 3].set_ylabel(L"$z$ (m)")
    ax[1, 3].set_title("buoyancy gradient")

    ax[2, 1].set_xlabel(L"$u$ (m s$^{-1}$)")
    ax[2, 1].set_ylabel(L"$z$ (m)")
    ax[2, 1].set_title("cross-ridge velocity")

    ax[2, 2].set_xlabel(L"$v$ (m s$^{-1}$)")
    ax[2, 2].set_ylabel(L"$z$ (m)")
    ax[2, 2].set_title("along-ridge velocity")

    ax[2, 3].set_xlabel(L"$w$ (m s$^{-1}$)")
    ax[2, 3].set_ylabel(L"$z$ (m)")
    ax[2, 3].set_title("vertical velocity")

    subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.9, wspace=0.3, hspace=0.6)

    for a in ax
        a.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
    end
    axins21.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))

    # color map
    colors = pl.cm.viridis(range(1, 0, length=5))

    # zoomed z
    ax[2, 1].set_ylim([z[ix, 1], z[ix, 1] + 200])

    # plot data from `datafiles`
    for i=1:size(datafiles, 1)
        # load
        b, chi, û, v, U, t, L, H0, Pr, f, N, symmetry, κ = loadCheckpointRot(datafiles[i])

        # convert to physical coordinates 
        u, w = rotate(û)

        # gradient
        forcing = -ẑDerivative(b).*sinθ

        # stratification
        Bz = N^2*cosθ + ẑDerivative(b)

        # colors and labels
        label = string("Day ", Int64(round(t/86400)))
        c = colors[i, :]

        # plot
        ax[1, 1].plot(Bz[ix, :],  z[ix, :], c=c)
        ax[1, 2].plot(chi[ix, :], z[ix, :], c=c)
        ax[1, 2].axvline(U[ix], c=c, lw=1.0, ls="--")
        ax[1, 3].plot(forcing[ix, :],  z[ix, :], c=c)
        ax[2, 1].plot(u[ix, :],   z[ix, :], c=c)
        ax[2, 2].plot(v[ix, :],   z[ix, :], c=c)
        ax[2, 3].plot(w[ix, :],   z[ix, :], c=c, label=label)
        axins21.plot(u[ix, :],    z[ix, :], c=c)
    end

    ax[2, 3].legend()

    savefig("profiles.png", bbox="inches")
end

"""
    advectionPlot(datafiles, iξ)

Plot advection terms from HDF5 snapshot files of buoyancy in the `datafiles` list
at for ξ = ξ[iξ].
"""
function advectionProfilePlot(datafiles, iξ)
    fig, ax = subplots(1, 3, figsize=(3.404*3, 3.404/1.62), sharey=true)
    colors = pl.cm.viridis(range(1, 0, length=5))
    for i=1:size(datafiles, 1)
        # load
        b, chi, û, v, U, t, L, H0, Pr, f, N, symmetry, κ = loadCheckpointRot(datafiles[i])
    
        adv = -N^2*û.*sinθ
        diff = ẑDerivative(κ.*(N^2*cosθ .+ ẑDerivative(b)))
        sum = adv + diff
    
        ax[1].plot(adv[iξ, :],  ẑ[iξ, :], c=colors[i, :], label=string("Day ", Int64(t/86400)))
        ax[2].plot(diff[iξ, :], ẑ[iξ, :], c=colors[i, :])
        ax[3].plot(sum[iξ, :],  ẑ[iξ, :], c=colors[i, :])
        ax[1].set_xlabel(L"-\hat{u} N^2 \sin\theta")
        ax[2].set_xlabel(L"[\kappa(N^2\cos\theta + b_\hat{z})]_\hat{z}")
        ax[3].set_xlabel(L"b_t")
    end
    ax[1].set_ylabel(L"\hat{z}")
    ax[1].legend()
    tight_layout()
    savefig("evolProfiles.png")
end

"""
    plotCurrentState(t, chi, û, v, b, iImg)

Generate ridgePlot's for each variable.
"""
function plotCurrentState(t, chi, û, v, b, iImg)
    # rotate 
    u, w = rotate(û)

    # plots
    ridgePlot(chi, b, @sprintf("streamfunction at t = %4d days", t/86400), L"$\chi$ (m$^2$ s$^{-1}$)")
    savefig(@sprintf("chi%03d.png", iImg))
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
