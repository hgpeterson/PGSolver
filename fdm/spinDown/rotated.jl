################################################################################
# Utility functions for rotated coordinates
################################################################################

"""
    u, w = rotate(û)

Rotate `û` into physical coordinate components `u` and `w`.
"""
function rotate(û)
    u = û*cos(θ)
    w = û*sin(θ)
    return u, w
end

"""
    saveCheckpointSpinDown(û, v, b, Px, t, i)

Save .h5 checkpoint file for state.
"""
function saveCheckpointSpinDown(û, v, b, Px, t, i)
    savefile = @sprintf("checkpoint%d.h5", i)
    file = h5open(savefile, "w")
    write(file, "û", û)
    write(file, "v", v)
    write(file, "b", b)
    write(file, "Px", Px)
    write(file, "t", t)
    write(file, "H0", H0)
    write(file, "Pr", Pr)
    write(file, "f", f)
    write(file, "N", N)
    write(file, "θ", θ)
    write(file, "canonical", canonical)
    write(file, "bottomIntense", bottomIntense)
    write(file, "κ", κ)
    write(file, "κ0", κ0)
    write(file, "κ1", κ1)
    write(file, "h", h)
    write(file, "α", α)
    close(file)
    println(savefile)
end

"""
    û, v, b, Px, t, L, H0, Pr, f, N, θ, canonical, bottomIntense, κ, κ0, κ1, h, α
        = loadCheckpointSpinDown(filename)

Load .h5 checkpoint file given by `filename`.
"""
function loadCheckpointSpinDown(filename)
    file = h5open(filename, "r")
    û = read(file, "û")
    v = read(file, "v")
    b = read(file, "b")
    Px = read(file, "Px")
    t = read(file, "t")
    H0 = read(file, "H0")
    Pr = read(file, "Pr")
    f = read(file, "f")
    N = read(file, "N")
    θ = read(file, "θ")
    canonical = read(file, "canonical")
    bottomIntense = read(file, "bottomIntense")
    κ = read(file, "κ")
    κ0 = read(file, "κ0")
    κ1 = read(file, "κ1")
    h = read(file, "h")
    α = read(file, "α")
    close(file)
    return û, v, b, Px, t, L, H0, Pr, f, N, θ, canonical, bottomIntense, κ, κ0, κ1, h, α
end
