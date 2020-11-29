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

        # canonical solution
        b, chi, û, v, U, t, L, H0, Pr, f, N, symmetry, κ = loadCheckpoint(string(folder, "canonical/checkpoint", tDay, ".h5"))
        ax.plot(v[ix, :], z[ix, :], label="canonical")

        # fixed solution
        b, chi, û, v, U, t, L, H0, Pr, f, N, symmetry, κ = loadCheckpoint(string(folder, "fixed/checkpoint", tDay, ".h5"))
        ax.plot(v[ix, :], z[ix, :], label="fixed")

        ax.legend(loc="upper right")
        tight_layout()
        fname = @sprintf("vProfiles%05d.png", tDay)
        savefig(fname)
        println(fname)
        close()
    end
end

function uProfile(folder)
    ix = 1
    tDay = 1000
    
    # setup plot
    fig, ax = subplots(1)
    ax.set_xlabel(L"cross-slope velocity, $u$ (m s$^{-1}$)")
    ax.set_ylabel(L"$z$ (m)")
    #= ax.set_xlim([-0.02, 0.02]) =#
    ax.set_ylim([z[ix, 1], z[ix, 1] + 100])
    ax.set_title(string(L"$t = $", tDay, " days"))

    # canonical solution
    b, chi, û, v, U, t, L, H0, Pr, f, N, symmetry, κ = loadCheckpoint(string(folder, "canonical/checkpoint", tDay, ".h5"))
    u, w = rotate(û)
    ax.plot(u[ix, :], z[ix, :], label="canonical")

    # fixed solution
    b, chi, û, v, U, t, L, H0, Pr, f, N, symmetry, κ = loadCheckpoint(string(folder, "fixed/checkpoint", tDay, ".h5"))
    u, w = rotate(û)
    ax.plot(u[ix, :], z[ix, :], label="fixed")

    ax.legend(loc="upper right")
    tight_layout()
    fname = @sprintf("uProfiles%05d.png", tDay)
    savefig(fname)
    println(fname)
    close()
end
