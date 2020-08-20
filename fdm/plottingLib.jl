################################################################################
# Functions useful for plotting
################################################################################

using PyPlot

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
    img = ax.pcolormesh(x, z, field, cmap="RdBu_r", vmin=vmin, vmax=vmax)
    colorbar(img, ax=ax, label=cbarLabel)

    # isopycnal contours
    ax.contour(x, z, B, 20, colors="k", alpha=0.3, linestyles="-")

    # ridge shading
    ax.fill_between(x[:, 1], z[:, 1], minimum(z), color="k", alpha=0.3)

    # labels
    ax.set_title(titleString)
    ax.set_xlabel(L"$x$ (m)")
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

function profilePlot(ax, u, v, w, b, iξ, day)
    bz = zDerivativeTF(b)

    ax[1, 1].plot(u[iξ, :], z[iξ, :], label=string("Day ", Int64(round(day))))
    ax[1, 2].plot(v[iξ, :], z[iξ, :])
    ax[2, 1].plot(w[iξ, :], z[iξ, :])
    ax[2, 2].plot(N^2 .+ bz[iξ, :], z[iξ, :])
end
