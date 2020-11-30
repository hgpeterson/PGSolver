################################################################################
# Utility functions for rotated coordinates
################################################################################

"""
    fẑ = ẑDerivative(field)

Compute dẑ(field) over 2D domain.
"""
function ẑDerivative(field)
    # allocate
    fẑ = zeros(nx, nz)

    # dẑ(field)
    for i=1:nx
        fẑ[i, :] .+= differentiate(field[i, :], ẑ[i, :])
    end

    return fẑ
end

"""
    u, w = rotate(û)

Rotate `û` into physical coordinate components `u` and `w`.
"""
function rotate(û)
    u = @. û*cosθ
    w = @. û*sinθ
    return u, w
end

"""
    saveCheckpointRot(b, chi, û, v, U, t)

Save .h5 checkpoint file for state `b` at time `t`.
"""
function saveCheckpointRot(b, chi, û, v, U, t)
    tDays = t/86400
    savefile = @sprintf("checkpoint%d.h5", tDays)
    file = h5open(savefile, "w")
    write(file, "b", b)
    write(file, "chi", chi)
    write(file, "û", û)
    write(file, "v", v)
    write(file, "U", U)
    write(file, "t", t)
    write(file, "L", L)
    write(file, "H0", H0)
    write(file, "Pr", Pr)
    write(file, "f", f)
    write(file, "N", N)
    write(file, "symmetry", symmetry)
    write(file, "κ", κ)
    close(file)
    println(savefile)
end

"""
    b, chi, û, v, U, t, L, H0, Pr, f, N, symmetry, κ = loadCheckpointRot(filename)

Load .h5 checkpoint file given by `filename`.
"""
function loadCheckpointRot(filename)
    file = h5open(filename, "r")
    b = read(file, "b")
    chi = read(file, "chi")
    û = read(file, "û")
    v = read(file, "v")
    U = read(file, "U")
    t = read(file, "t")
    L = read(file, "L")
    H0 = read(file, "H0")
    Pr = read(file, "Pr")
    f = read(file, "f")
    N = read(file, "N")
    symmetry = read(file, "symmetry")
    κ = read(file, "κ")
    close(file)
    return b, chi, û, v, U, t, L, H0, Pr, f, N, symmetry, κ
end
