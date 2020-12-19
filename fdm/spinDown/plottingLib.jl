################################################################################
# Functions useful for plotting
################################################################################

pl = pyimport("matplotlib.pylab")
inset_locator = pyimport("mpl_toolkits.axes_grid1.inset_locator")

"""
    profilePlot(datafiles)

Plot profiles from HDF5 snapshot files the `datafiles` list.
"""
function profilePlot(datafiles)
    # init plot
    fig, ax = subplots(1, 3, figsize=(3.404*3, 3.404/1.62))

    # insets
    axins21 = inset_locator.inset_axes(ax[2], width="40%", height="40%")

    ax[1].set_xlabel(L"$B_\hat{z}$ (s$^{-2}$)")
    ax[1].set_ylabel(L"$\hat{z}$ (m)")
    ax[1].set_title("stratification")

    ax[2].set_xlabel(L"$\hat{u}$ (m s$^{-1}$)")
    ax[2].set_ylabel(L"$\hat{z}$ (m)")
    ax[2].set_title("cross-ridge velocity")

    ax[3].set_xlabel(L"$v$ (m s$^{-1}$)")
    ax[3].set_ylabel(L"$\hat{z}$ (m)")
    ax[3].set_title("along-ridge velocity")

    subplots_adjust(bottom=0.2, top=0.9, left=0.1, right=0.9, wspace=0.3, hspace=0.6)

    for a in ax
        a.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
    end
    axins21.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))

    # color map
    colors = pl.cm.viridis(range(1, 0, length=size(datafiles, 1)-1))

    # zoomed z
    ax[1].set_ylim([ẑ[1], ẑ[1] + 100])
    ax[2].set_ylim([ẑ[1], ẑ[1] + 100])
    ax[3].set_ylim([ẑ[1], ẑ[1] + 100])

    # plot data from `datafiles`
    for i=1:size(datafiles, 1)
        # load
        û, v, b, Px, t, L, H0, Pr, f, N, θ, canonical, bottomIntense, κ, κ0, κ1, h, α = loadCheckpointSpinDown(datafiles[i])

        # stratification
        Bẑ = N^2*cos(θ) .+ differentiate(b, ẑ)

        # colors and labels
        label = string("Day ", Int64(round(t/86400)))
        if i == 1
            c = "k"
        else
            c = colors[i-1, :]
        end

        # plot
        ax[1].plot(Bẑ, ẑ, c=c, label=label)
        ax[2].plot(û,  ẑ, c=c)
        ax[3].plot(v,  ẑ, c=c)
        ax[3].axvline(Px/f, lw=1.0, c=c, ls="--")
        axins21.plot(û, ẑ, c=c)
    end

    ax[1].legend()

    savefig("profiles.png", bbox="inches")
    close()
end
