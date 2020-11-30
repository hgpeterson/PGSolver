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

function canonicalSteadyTheory(z, κ0, κ1, h, N, f, θ, Pr)
    S = N^2*tan(θ)^2/f^2
    q = (f^2*cos(θ)^2*(1 + S*Pr)/(4*Pr*(κ0 + κ1)^2))^(1/4)
    Bz = @. N^2*cos(θ)*(κ0/(κ0 + κ1*exp(-z/h)) + κ1*exp(-z/h)/(κ0 + κ1*exp(-z/h))*S*Pr/(1 + S*Pr) - (κ0/(κ0 + κ1) + κ1/(κ0 + κ1)*S*Pr/(1 + S*Pr))*exp(-q*z)*(cos(q*z) + sin(q*z)))
    u = @. -κ1*cot(θ)*exp(-z/h)/h*S*Pr/(1 + S*Pr) + 2*q*cot(θ)*(κ0 + κ1*S*Pr/(1 + S*Pr))*exp(-q*z)*sin(q*z)
    vz = @. 1 # FIXME
    v = cumtrapz(vz, z)
    return Bz, u, v
end
