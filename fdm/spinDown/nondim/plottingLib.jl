################################################################################
# Functions useful for plotting
################################################################################

pl = pyimport("matplotlib.pylab")
inset_locator = pyimport("mpl_toolkits.axes_grid1.inset_locator")

"""
    profilePlot(datafiles)

Plot profiles from HDF5 snapshot files the `datafiles` list.
"""
function profilePlot(datafiles; fname="profiles.png")
    # init plot
    fig, ax = subplots(1, 3, figsize=(3.404*3, 3.404/1.62), sharey=true)

    # insets
    #= axins21 = inset_locator.inset_axes(ax[2], width="40%", height="40%") =#

    ax[1].set_xlabel(L"$B_\hat{z}$")
    ax[1].set_ylabel(L"$\hat{z}$")
    ax[1].set_title("stratification")

    ax[2].set_xlabel(L"$\hat{u}$")
    ax[2].set_title("cross-ridge velocity")

    ax[3].set_xlabel(L"$v$")
    ax[3].set_title("along-ridge velocity")

    subplots_adjust(bottom=0.2, top=0.9, left=0.1, right=0.9, wspace=0.3, hspace=0.6)

    for a in ax
        a.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
    end
    #= axins21.ticklabel_format(style="sci", axis="x", scilimits=(0, 0)) =#

    # color map
    colors = pl.cm.viridis(range(1, 0, length=size(datafiles, 1)-1))

    # zoomed z
    ax[1].set_ylim([0, 10])

    # plot data from `datafiles`
    for i=1:size(datafiles, 1)
        # load
        û, v, b, Px, t, H, Pr, S, canonical, bottomIntense, κ, κ0, κ1, h, α = loadCheckpointSpinDown(datafiles[i])

        # stratification
        Bẑ = 1 .+ differentiate(b, ẑ)

        # colors and labels
        label = string(L"$t/\tau_A$ = ", Int64(round(t/τ_A)))
        if i == 1
            c = "k"
        else
            c = colors[i-1, :]
        end

        # plot
        ax[1].plot(Bẑ, ẑ, c=c, label=label)
        ax[2].plot(û,  ẑ, c=c)
        ax[3].plot(v,  ẑ, c=c)
        ax[3].axvline(Px, lw=1.0, c=c, ls="--")
        #= axins21.plot(û, ẑ, c=c) =#
    end

    ax[1].legend()

    savefig(fname)
    close()
end

"""
    MR91plot(datafiles)

Plot profiles from HDF5 snapshot files the `datafiles` list as in MR91.
"""
function MR91plot(datafiles)
    # init plot
    fig, ax = subplots(1, 2, figsize=(3.404*2, 3.404/1.62), sharey=true)

    ax[1].set_xlabel(L"$\hat{u}, \hat{v}$ (cm s$^{-1}$)")
    ax[1].set_ylabel(L"$\hat{z}$ (cm)")
    ax[1].set_title("velocities")

    ax[2].set_xlabel(L"$\overline{\rho} + \rho'$ (g cm$^{-3}$)")
    ax[2].set_title("density")

    subplots_adjust(bottom=0.2, top=0.9, left=0.1, right=0.9, wspace=0.3, hspace=0.6)

    for a in ax
        a.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
    end

    # color map
    #= colors = pl.cm.viridis(range(1, 0, length=size(datafiles, 1)-1)) =#
    colors = ["tab:blue", "tab:orange", "tab:green"]

    # xlim, ylim
    ax[1].set_xlim([-1.2, 0.4])
    ax[1].set_ylim([0, 2])
    ax[1].set_xticks(-1.2:0.2:0.4)
    ax[2].set_xlim([-0.01, 0.01])

    # plot data from `datafiles`
    for i=1:size(datafiles, 1)
        # load
        û, v, b, Px, t, L, H0, Pr, f, N, θ, canonical, bottomIntense, κ, κ0, κ1, h, α = loadCheckpointSpinDown(datafiles[i])

        # colors and labels
        #= label = string("Day ", Int64(round(t/86400))) =#
        label = string(Int64(round(t)), " s")
        c = colors[i]

        # plot
        ax[1].plot(100*û, 100*(ẑ .+ H0), c=c, ls="-")
        ax[1].plot(100*v, 100*(ẑ .+ H0), c=c, ls="--")
        ax[2].plot(-(b + N^2*(ẑ .+ H0))/9.81*1, 100*(ẑ .+ H0), c=c, label=label)
    end

    ax[2].legend()

    savefig("profiles.png", bbox="inches")
    close()
end

"""
    canonicalVsNoTransport(folder)

Plot profiles from HDF5 snapshot files in `folder` for both canonical and no_transport case.
"""
function canonicalVsNoTransport(folder)
    # init plot
    fig, ax = subplots(3, 2, figsize=(3.404*2, 3*3.404/1.62), sharey=true)

    ax[1, 1].set_xlabel(L"cross-ridge flow, $\hat{u}$ (m s$^{-1}$)")
    ax[1, 1].set_ylabel(L"$\hat{z}$ (m)")
    ax[1, 1].set_title("canonical")

    ax[1, 2].set_xlabel(L"cross-ridge flow, $\hat{u}$ (m s$^{-1}$)")
    ax[1, 2].set_title("transport-constrained")

    ax[2, 1].set_xlabel(L"along-ridge flow, $\hat{v}$ (m s$^{-1}$)")
    ax[2, 1].set_ylabel(L"$\hat{z}$ (m)")

    ax[2, 2].set_xlabel(L"along-ridge flow, $\hat{v}$ (m s$^{-1}$)")

    ax[3, 1].set_xlabel(L"stratification, $B_\hat{z}$ (s$^{-2}$)")
    ax[3, 1].set_ylabel(L"$\hat{z}$ (m)")

    ax[3, 2].set_xlabel(L"stratification, $B_\hat{z}$ (s$^{-2}$)")

    for a in ax
        a.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
    end

    # color map
    colors = pl.cm.viridis(range(1, 0, length=5))

    # zoomed z
    ax[1, 1].set_ylim([0, 100])
    ax[2, 1].set_ylim([0, 100])
    ax[3, 1].set_ylim([0, 100])

    # plot data from folder
    cases = ["canonical", "no_transport"]
    S = 0.0 # need S later for plot
    for j=1:2
        case = cases[j]
        for i=0:5
            # load
            û, v, b, Px, t, L, H0, Pr, f, N, θ, canonical, bottomIntense, κ, κ0, κ1, h, α = loadCheckpointSpinDown(string(folder, "/", case, "/checkpoint", i, ".h5"))
            S = N^2*tan(θ)^2/f^2
            τ = 1/abs(S*f*cos(θ))

            # stratification
            Bẑ = N^2*cos(θ) .+ differentiate(b, ẑ)

            # colors and labels
            if i == 0
                c = "k"
                label = L"$t = 0$"
            else
                c = colors[i, :]
                label = string(L"$t/\tau$ = ", Int64(round(t/τ)))
            end

            # plot
            ax[1, j].plot(û,  ẑ, c=c, label=label)
            ax[2, j].plot(v,  ẑ, c=c, label=label)
            ax[3, j].plot(Bẑ, ẑ, c=c, label=label)
            ax[2, j].axvline(Px/f, lw=1.0, c=c, ls="--")
        end
    end

    ax[1, 1].legend()

    ax[1, 1].annotate("(a)", (0.06, 0.92), xycoords="axes fraction")
    ax[1, 2].annotate("(b)", (0.06, 0.92), xycoords="axes fraction")
    ax[2, 1].annotate("(c)", (0.06, 0.92), xycoords="axes fraction")
    ax[2, 2].annotate("(d)", (0.06, 0.92), xycoords="axes fraction")
    ax[3, 1].annotate("(e)", (0.06, 0.92), xycoords="axes fraction")
    ax[3, 2].annotate("(f)", (0.06, 0.92), xycoords="axes fraction")

    ax[1, 2].annotate(string(L"$S \approx $", @sprintf("%1.1e", S)),  (0.6, 0.6), xycoords="axes fraction", size=10)
    ax[1, 2].annotate(string(L"$H = $", @sprintf("%1d km", H0/1000)), (0.6, 0.5), xycoords="axes fraction", size=10)

    #= ax[2, 2].annotate(L"$P_x/f$", xy=(0.18, 0.5), xytext=(0.03, 0.6), xycoords="axes fraction", arrowprops=Dict("arrowstyle" => "->")) =#

    tight_layout()
    savefig("fig.pdf")
    close()
end
