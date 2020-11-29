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
