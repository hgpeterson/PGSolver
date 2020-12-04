# for loading rotated data
include("rotatedCoords/rotated.jl")

function vAnimation(folder)
    ix = 1
    tDays = 0:10:1000
    for tDay in tDays
        # setup plot
        fig, ax = subplots(1)
        ax.set_xlabel(L"along-slope velocity, $v$ (m s$^{-1}$)")
        ax.set_ylabel(L"$z$ (m)")
        #= ax.set_xlim([-0.02, 0.02]) =#
        ax.set_xlim([-0.03, 0.03])
        ax.set_title(string(L"$t = $", tDay, " days"))

        # full 2D solution
        b, chi, uξ, uη, uσ, U, t, L, H0, Pr, f, N, symmetry, ξVariation, κ = loadCheckpointTF(string(folder, "full2D/checkpoint", tDay, ".h5"))
        u, v, w = transformFromTF(uξ, uη, uσ)
        ax.plot(v[ix, :], z[ix, :], label="full 2D")

        # canonical 1D solution
        b, chi, û, v, U, t, L, H0, Pr, f, N, symmetry, κ = loadCheckpointRot(string(folder, "canonical1D/checkpoint", tDay, ".h5"))
        ax.plot(v[ix, :], z[ix, :], label="canonical 1D")

        # fixed 1D solution
        b, chi, û, v, U, t, L, H0, Pr, f, N, symmetry, κ = loadCheckpointRot(string(folder, "fixed1D/checkpoint", tDay, ".h5"))
        ax.plot(v[ix, :], z[ix, :], "k--", label="fixed 1D")

        ax.legend(loc="upper right")
        tight_layout()
        fname = @sprintf("vProfiles%05d.png", tDay)
        savefig(fname)
        println(fname)
        close()
    end
end

function uAnimation(folder)
    ix = 1
    tDays = 0:10:1000
    #= tDays = 1000 =#
    for tDay in tDays
        # setup plot
        fig, ax = subplots(1)
        ax.set_xlabel(L"cross-slope velocity, $u$ (m s$^{-1}$)")
        ax.set_ylabel(L"$z$ (m)")
        ax.set_xlim([-0.0005, 0.007])
        ax.set_ylim([z[ix, 1], z[ix, 1] + 200])
        ax.set_title(string(L"$t = $", tDay, " days"))

        # full 2D solution
        b, chi, uξ, uη, uσ, U, t, L, H0, Pr, f, N, symmetry, ξVariation, κ = loadCheckpointTF(string(folder, "full2D/checkpoint", tDay, ".h5"))
        u, v, w = transformFromTF(uξ, uη, uσ)
        ax.plot(u[ix, :], z[ix, :], label="full 2D")

        # canonical 1D solution
        b, chi, û, v, U, t, L, H0, Pr, f, N, symmetry, κ = loadCheckpointRot(string(folder, "canonical1D/checkpoint", tDay, ".h5"))
        u, w = rotate(û)
        ax.plot(u[ix, :], z[ix, :], label="canonical 1D")
        #= # steady state =#
        #= if !bottomIntense =#
        #=     Bz, u, v = canonicalSteadyTheory(z[ix, :] .+ H0, κ0, κ1, h, N, f, θ[ix], Pr) =#
        #= else =#
        #=     Bz, u, v = canonicalSteadyTheory(z[ix, :] .+ H0, κ1, 0, 1, N, f, θ[ix], Pr) =#
        #= end =#
        #= ax.plot(u, z[ix, :], c="tab:green", ls="--", label="canonical 1D steady state") =#
        

        # fixed 1D solution
        b, chi, û, v, U, t, L, H0, Pr, f, N, symmetry, κ = loadCheckpointRot(string(folder, "fixed1D/checkpoint", tDay, ".h5"))
        u, w = rotate(û)
        ax.plot(u[ix, :], z[ix, :], "k--", label="fixed 1D")

        ax.legend(loc="upper right")
        tight_layout()
        fname = @sprintf("uProfiles%05d.png", tDay)
        savefig(fname)
        println(fname)
        close()
    end
end

function uvAnimation(folder)
    ix = 1
    tDays = 0:10:1000
    #= tDays = 1000 =#
    for tDay in tDays
        # setup plot
        fig, ax = subplots(1, 2, figsize=(6.5, 6.5/1.62/2))

        ax[1].set_xlabel(L"cross-slope velocity, $u$ (m s$^{-1}$)")
        ax[1].set_ylabel(L"$z$ (m)")
        ax[1].set_xlim([-0.0005, 0.007])
        ax[1].set_ylim([z[ix, 1], z[ix, 1] + 200])
        ax[1].set_title(string(L"$t = $", tDay, " days"))

        ax[2].set_xlabel(L"along-slope velocity, $v$ (m s$^{-1}$)")
        ax[2].set_ylabel(L"$z$ (m)")
        ax[2].set_xlim([-0.03, 0.03])
        ax[2].set_title(string(L"$t = $", tDay, " days"))

        # full 2D solution
        b, chi, uξ, uη, uσ, U, t, L, H0, Pr, f, N, symmetry, ξVariation, κ = loadCheckpointTF(string(folder, "full2D/checkpoint", tDay, ".h5"))
        u, v, w = transformFromTF(uξ, uη, uσ)
        ax[1].plot(u[ix, :], z[ix, :], label="full 2D")
        ax[2].plot(v[ix, :], z[ix, :], label="full 2D")

        # canonical 1D solution
        b, chi, û, v, U, t, L, H0, Pr, f, N, symmetry, κ = loadCheckpointRot(string(folder, "canonical1D/checkpoint", tDay, ".h5"))
        u, w = rotate(û)
        ax[1].plot(u[ix, :], z[ix, :], label="canonical 1D")
        ax[2].plot(v[ix, :], z[ix, :], label="canonical 1D")

        # fixed 1D solution
        b, chi, û, v, U, t, L, H0, Pr, f, N, symmetry, κ = loadCheckpointRot(string(folder, "fixed1D/checkpoint", tDay, ".h5"))
        u, w = rotate(û)
        ax[1].plot(u[ix, :], z[ix, :], "k--", label="fixed 1D")
        ax[2].plot(v[ix, :], z[ix, :], "k--", label="fixed 1D")

        ax[1].legend(loc="upper right")
        tight_layout()
        fname = @sprintf("uvProfiles%04d.png", tDay)
        #= fname = @sprintf("uvProfiles%04d.pdf", tDay) =#
        savefig(fname, dpi=200)
        println(fname)
        close()
    end
end

function canonicalSteadyTheory(z, κ0, κ1, h, N, f, θ, Pr)
    S = N^2*tan(θ)^2/f^2
    q = (f^2*cos(θ)^2*(1 + S*Pr)/(4*Pr*(κ0 + κ1)^2))^(1/4)
    Bz = @. N^2*cos(θ)*(κ0/(κ0 + κ1*exp(-z/h)) + κ1*exp(-z/h)/(κ0 + κ1*exp(-z/h))*S*Pr/(1 + S*Pr) - (κ0/(κ0 + κ1) + κ1/(κ0 + κ1)*S*Pr/(1 + S*Pr))*exp(-q*z)*(cos(q*z) + sin(q*z)))
    u = @. -κ1*cot(θ)*exp(-z/h)/h*S*Pr/(1 + S*Pr) + 2*q*cot(θ)*(κ0 + κ1*S*Pr/(1 + S*Pr))*exp(-q*z)*sin(q*z)
    vz = @. 1 # FIXME
    v = cumtrapz(vz, z)
    return Bz, u, v
end

function idealRidge()
    fig, ax = subplots(1)
    ax.fill_between(x[:, 1]/1000, z[:, 1], minimum(z), color="k", alpha=0.3, lw=0.0)
    ax.annotate("", xy=(1500, -2800), xytext=(1500, 0), xycoords="data", arrowprops=Dict("arrowstyle" => "<->"))
    ax.annotate(L"$H(x)$", (1600, -1000), xycoords="data")
    ax.annotate(L"$κ(z)$", (1000, -800), xycoords="data")
    ax.annotate(L"$\bf{u}$", (1200, -2200), xycoords="data")
    ax.annotate("periodic", xy=(0, -400), xytext=(300, -450), xycoords="data", arrowprops=Dict("arrowstyle" => "<->"))
    ax.spines["left"].set_visible(false)
    ax.spines["bottom"].set_visible(false)
    ax.set_ylim([minimum(z), 0])
    ax.set_xlabel(L"$x$ (km)")
    ax.set_ylabel(L"$z$ (m)")
    tight_layout()
    savefig("ideal_ridge.svg")
end

function uBalance(folder)
    iξ = 1

    fig, ax = subplots(1)

    # full 2D
    b, chi, uξ, uη, uσ, U, t, L, H0, Pr, f, N, symmetry, ξVariation, κ = loadCheckpointTF(string(folder, "full2D/checkpoint1000.h5"))
    u, v, w = transformFromTF(uξ, uη, uσ)
    û = @. cosθ*u + sinθ*w
    ẑ = @. z*cosθ
    ax.plot(û[iξ, :]/maximum(û[iξ, :]), ẑ[iξ, :], label="full 2D")

    b, chi, û, v, U, t, L, H0, Pr, f, N, symmetry, κ = loadCheckpointRot(string(folder, "canonical1D/checkpoint1000.h5"))
    ax.plot(û[iξ, :]/maximum(û[iξ, :]), ẑ[iξ, :], label="canonical 1D")

    ax.set_xlabel(L"$\hat{u}$ (rescaled)")
    ax.set_ylabel(L"$\hat{z}$")
    ax.set_ylim([ẑ[iξ, 1], ẑ[iξ, 1] + 500])
    ax.set_xticks([0])
    ax.set_yticks([ẑ[iξ, 1], ẑ[iξ, 1] + 45])
    ax.set_yticklabels([0, L"5$\delta$"])
    ax.legend(loc="upper right")
    ax.annotate("far-field downwelling", xy=(0.07, 0.4), xytext=(0.3, 0.5), xycoords="axes fraction", arrowprops=Dict("arrowstyle" => "->"))
    tight_layout()
    savefig("ubalance.pdf")
end

function uvRidge(folder)
    b, chi, uξ, uη, uσ, U, t, L, H0, Pr, f, N, symmetry, ξVariation, κ = loadCheckpointTF(string(folder, "full2D/checkpoint1000.h5"))
    u, v, w = transformFromTF(uξ, uη, uσ)

    ridgePlot(u, b, "cross-ridge velocity", L"$u$ (m s$^{-1}$)"; vext=1e-5)
    savefig("u1000.pdf")
    ridgePlot(v, b, "along-ridge velocity", L"$v$ (m s$^{-1}$)")
    savefig("v1000.pdf")
    close("all")
end

function uvPrScaling(folder)
    ix = 1
    tDay = 1000
    Prs = [1, 10, 100, 1000]
    #= lss = ["-", "--", "-.", ":"] =#
    alphas = [1.0, 0.75, 0.5, 0.25]

    # setup plot
    fig, ax = subplots(1, 2, figsize=(6.5, 6.5/1.62/2))
    ax[1].set_xlabel(L"cross-slope velocity, $u$ (m s$^{-1}$)")
    ax[1].set_ylabel(L"$z$ (m)")
    ax[1].set_xlim([-0.0005, 0.0035])
    ax[1].set_ylim([z[ix, 1], z[ix, 1] + 200])
    ax[1].set_title(string(L"$t = $", tDay, " days"))
    ax[2].set_xlabel(L"along-slope velocity, $v$ (m s$^{-1}$)")
    ax[2].set_ylabel(L"$z$ (m)")
    ax[2].set_xlim([-0.015, 0.015])
    ax[2].set_title(string(L"$t = $", tDay, " days"))

    # loop
    for i=1:size(Prs, 1)
        σ = Prs[i]
        ls = "-"
        alpha = alphas[i]

        # full 2D solution
        b, chi, uξ, uη, uσ, U, t, L, H0, Pr, f, N, symmetry, ξVariation, κ = loadCheckpointTF(string(folder, "full2D/Pr", σ, "/checkpoint", tDay, ".h5"))
        u, v, w = transformFromTF(uξ, uη, uσ)
        ax[1].plot(u[ix, :], z[ix, :], c="tab:blue", ls=ls, label="full 2D", alpha=alpha)
        ax[2].plot(v[ix, :], z[ix, :], c="tab:blue", ls=ls, label="full 2D", alpha=alpha)

        # canonical 1D solution
        b, chi, û, v, U, t, L, H0, Pr, f, N, symmetry, κ = loadCheckpointRot(string(folder, "canonical1D/Pr", σ, "/checkpoint", tDay, ".h5"))
        u, w = rotate(û)
        ax[1].plot(u[ix, :], z[ix, :], c="tab:orange", ls=ls, label="canonical 1D", alpha=alpha)
        ax[2].plot(v[ix, :], z[ix, :], c="tab:orange", ls=ls, label="canonical 1D", alpha=alpha)

        # fixed 1D solution
        b, chi, û, v, U, t, L, H0, Pr, f, N, symmetry, κ = loadCheckpointRot(string(folder, "fixed1D/Pr", σ, "/checkpoint", tDay, ".h5"))
        u, w = rotate(û)
        ax[1].plot(u[ix, :], z[ix, :], c="k", ls="--", label="fixed 1D", alpha=alpha)
        ax[2].plot(v[ix, :], z[ix, :], c="k", ls="--", label="fixed 1D", alpha=alpha)

        if i == 1
            ax[1].legend(loc="upper right")
        end
    end
    tight_layout()
    fname = "uvPrScaling.pdf"
    savefig(fname)
    println(fname)
    close()
end

function BzChiPrScaling(folder)
    ix = 1
    tDay = 1000
    Prs = [1, 10, 100, 1000]
    #= lss = ["-", "--", "-.", ":"] =#
    alphas = [1.0, 0.75, 0.5, 0.25]

    # setup plot
    fig, ax = subplots(1, 2, figsize=(6.5, 6.5/1.62/2))
    ax[1].set_xlabel(L"stratification, $B_z$ (s$^{-2}$)")
    ax[1].set_ylabel(L"$z$ (m)")
    #= ax[1].set_xlim([-0.0005, 0.0035]) =#
    ax[1].set_title(string(L"$t = $", tDay, " days"))
    ax[1].ticklabel_format(style="sci", scilimits=(-4, 4))
    ax[2].set_xlabel(L"streamfunction, $\chi$ (m$^2$ s$^{-1}$)")
    ax[2].set_ylabel(L"$z$ (m)")
    ax[2].set_xlim([-0.02, 0.18])
    ax[2].set_title(string(L"$t = $", tDay, " days"))

    # loop
    for i=1:size(Prs, 1)
        σ = Prs[i]
        ls = "-"
        alpha = alphas[i]

        # full 2D solution
        b, chi, uξ, uη, uσ, U, t, L, H0, Pr, f, N, symmetry, ξVariation, κ = loadCheckpointTF(string(folder, "full2D/Pr", σ, "/checkpoint", tDay, ".h5"))
        Bz = N^2 .+ zDerivativeTF(b)
        ax[1].plot(Bz[ix, :], z[ix, :], c="tab:blue", ls=ls, label="full 2D", alpha=alpha)
        ax[2].plot(chi[ix, :], z[ix, :], c="tab:blue", ls=ls, label="full 2D", alpha=alpha)

        # canonical 1D solution
        b, chi, û, v, U, t, L, H0, Pr, f, N, symmetry, κ = loadCheckpointRot(string(folder, "canonical1D/Pr", σ, "/checkpoint", tDay, ".h5"))
        Bz = N^2*cosθ[ix, :] .+ differentiate(b[ix, :], z[ix, :].*cosθ[ix, :])
        ax[1].plot(Bz, z[ix, :], c="tab:orange", ls=ls, label="canonical 1D", alpha=alpha)
        ax[2].plot(chi[ix, :], z[ix, :], c="tab:orange", ls=ls, label="canonical 1D", alpha=alpha)

        # fixed 1D solution
        b, chi, û, v, U, t, L, H0, Pr, f, N, symmetry, κ = loadCheckpointRot(string(folder, "fixed1D/Pr", σ, "/checkpoint", tDay, ".h5"))
        Bz = N^2*cosθ[ix, :] .+ differentiate(b[ix, :], z[ix, :].*cosθ[ix, :])
        ax[1].plot(Bz, z[ix, :], c="k", ls="--", label="fixed 1D", alpha=alpha)
        ax[2].plot(chi[ix, :], z[ix, :], c="k", ls="--", label="fixed 1D", alpha=alpha)

        if i == 1
            ax[1].legend(loc="upper left")
        end
    end
    tight_layout()
    fname = "BzChiPrScaling.pdf"
    savefig(fname)
    println(fname)
    close()
end

function pressureRidgePlots(dfile)
    # read
    b, chi, uξ, uη, uσ, U, t, L, H0, Pr, f, N, symmetry, ξVariation, κ = loadCheckpointTF(dfile)

    # compute p_x
    u, v, w = transformFromTF(uξ, uη, uσ)
    px = f*v + zDerivativeTF(Pr*κ.*zDerivativeTF(u)) 

    # compute p
    p = zeros(size(b))
    p[:, end] = cumtrapz(px[:, end], ξ) # assume p = 0 at top left
    for i=1:nξ
        hydrostatic = H(ξ[i])*cumtrapz(b[i, :], σ)
        p[i, :] = hydrostatic .+ (p[i, end] - hydrostatic[end]) # integration constant from int(px)
    end

    ridgePlot(p, b, "pressure", L"$p$ (m$^2$ s$^{-2}$)"; cmap="viridis")
    savefig("p1000.pdf")
    close()
    ridgePlot(px, b, "pressure gradient", L"$p_x$ (m s$^{-2}$)")
    savefig("px1000.pdf")
    close()
end
