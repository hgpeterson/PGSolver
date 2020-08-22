################################################################################
# Functions useful for plotting
################################################################################

using PyPlot, Printf

"""
    ax = ridgePlot(field, b, titleString, cbarLabel; vext)

Create 2D plot of `field` with isopycnals given by the buoyancy perturbation `b`.
Set the title to `titleString` and colorbar label to `cbarLabel`. Return the axis 
handle `ax`.

Optional: set the vmin/vmax manually with vext.
"""
function ridgePlot(field, b, titleString, cbarLabel; vext=nothing)
    # full buoyancy for isopycnals
    B = N^2*z + b 

    fig, ax = subplots(1)

    # set min and max
    if vext == nothing
        vmax = maximum(abs.(field))
        vmin = -vmax
    else
        vmax = vext
        vmin = -vext
    end

    # 2D plot
    img = ax.pcolormesh(r, z, field, cmap="RdBu_r", vmin=vmin, vmax=vmax)
    colorbar(img, ax=ax, label=cbarLabel)

    # isopycnal contours
    ax.contour(r, z, B, 20, colors="k", alpha=0.3, linestyles="-")

    # ridge shading
    ax.fill_between(r[:, 1], z[:, 1], minimum(z), color="k", alpha=0.3)

    # labels
    ax.set_title(titleString)
    ax.set_xlabel(L"$r$ (m)")
    ax.set_ylabel(L"$z$ (m)")
    
    return ax
end

function profilePlotInit()
    fig, ax = subplots(2, 2, sharey=true)

    ax[1, 1].set_xlabel(L"$u$, (m s$^{-1}$)")
    ax[1, 1].set_ylabel(L"$z$, (m)")
    ax[1, 1].set_title("cross-ridge velocity")

    ax[1, 2].set_xlabel(L"$v$, (m s$^{-1}$)")
    #= ax[1, 2].set_ylabel(L"$z$, (m)") =#
    ax[1, 2].set_title("along-ridge velocity")

    ax[2, 1].set_xlabel(L"$w$, (m s$^{-1}$)")
    ax[2, 1].set_ylabel(L"$z$, (m)")
    ax[2, 1].set_title("vertical velocity")

    ax[2, 2].set_xlabel(L"$B_z$, (s$^{-2}$)")
    #= ax[2, 2].set_ylabel(L"$z$, (m)") =#
    ax[2, 2].set_title("stratification")

    tight_layout()

    ax[1, 1].ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
    ax[1, 2].ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
    ax[2, 1].ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
    ax[2, 2].ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
    
    return ax
end

function profilePlot(ax, u, v, w, b, iρ, day)
    bz = zDerivativeTF(b)

    ax[1, 1].plot(u[iρ, :], z[iρ, :], label=string("Day ", Int64(round(day))))
    ax[1, 2].plot(v[iρ, :], z[iρ, :])
    ax[2, 1].plot(w[iρ, :], z[iρ, :])
    ax[2, 2].plot(N^2 .+ bz[iρ, :], z[iρ, :])
end

"""
    plotCurrentState(t, chi, chiEkman, u, v, w, b, iImg)

Plot the buoyancy and velocity state of the model at time `t` using label number `iImg`.
"""
function plotCurrentState(t, chi, chiEkman, u, v, w, b, iImg)
    ridgePlot(chi, b, @sprintf("streamfunction at t = %.1e s", t), L"$\chi$ (m$^2$ s$^{-1}$)")
    savefig(@sprintf("chi%03d.png", iImg))
    close()

    ridgePlot(chiEkman, b, @sprintf("streamfunction theory at t = %.1e s", t), L"$\chi$ (m$^2$ s$^{-1}$)")
    savefig(@sprintf("chiEkman%03d.png", iImg))
    close()

    ridgePlot(b, b, @sprintf("buoyancy perturbation at t = %.1e s", t), L"$b$ (m s$^{-2}$)")
    savefig(@sprintf("b%03d.png", iImg))
    close()

    #= ridgePlot(u, b, @sprintf("cross-ridge velocity at t = %.1f days", t/86400), L"$u$ (m s$^{-1}$)"; vext=5e-5) =#
    ridgePlot(u, b, @sprintf("cross-ridge velocity at t = %.1e s", t), L"$u$ (m s$^{-1}$)")
    savefig(@sprintf("u%03d.png", iImg))
    close()

    ridgePlot(v, b, @sprintf("along-ridge velocity at t = %.1e s", t), L"$v$ (m s$^{-1}$)")
    savefig(@sprintf("v%03d.png", iImg))
    close()

    ridgePlot(w, b, @sprintf("vertical velocity at t = %.1e s", t), L"$w$ (m s$^{-1}$)")
    savefig(@sprintf("w%03d.png", iImg))
    close()
end
