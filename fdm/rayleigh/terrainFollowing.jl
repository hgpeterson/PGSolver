"""
    fξ = ξDerivativeTF(field)

Compute dξ(`field`) in terrian-following coordinates.
"""
function ξDerivativeTF(field)
    # allocate
    fξ = zeros(nξ, nσ)

    # dξ(field)
    for j=1:nσ
        # use the fact that ξ is evenly spaced and periodic
        fξ[2:end-1, j] = (field[3:end, j] - field[1:end-2, j])/(2*dξ)
        fξ[1, j] = (field[2, j] - field[nξ, j])/(2*dξ)
        fξ[end, j] = (field[1, j] - field[end-1, j])/(2*dξ)
    end

    return fξ
end

"""
    fσ = σDerivativeTF(field)

Compute dσ(`field`) in terrian-following coordinates.
"""
function σDerivativeTF(field)
    # allocate
    fσ = zeros(nξ, nσ)

    # dσ(field)
    for i=1:nξ
        fσ[i, :] .+= differentiate(field[i, :], σ)
    end

    return fσ
end

"""
    fx = xDerivativeTF(field)

Compute dx(`field`) in terrian-following coordinates.
Note: dx() = dξ() - dx(H)*σ*dσ()/H
"""
function xDerivativeTF(field)
    # dξ(field)
    fx = ξDerivativeTF(field)

    # -dx(H)*σ*dσ(field)/H
    fx -= Hx.(x).*σσ.*σDerivativeTF(field)./H.(x)

    return fx
end

"""
    fz = zDerivativeTF(field)

Compute dz(`field`) in terrian-following coordinates.
Note: dz() = dσ()/H
"""
function zDerivativeTF(field)
    # dσ(field)/H
    fz = σDerivativeTF(field)./H.(x)

    return fz
end

"""
    u, v, w = transformFromTF(uξ, uη, uσ)

Transform from terrain-following coordinates to cartesian coordinates.
"""
function transformFromTF(uξ, uη, uσ)
    u = uξ

    v = uη

    w = uσ.*H.(x) + σσ.*Hx.(x).*u

    return u, v, w
end

"""
    uξ, uη, uσ = transformToTF(u, v, w)

Transform from cartesian coordinates to terrain-following coordinates.
"""
function transformToTF(u, v, w)
    uξ = u

    uη = v

    uσ = (w - σσ.*Hx.(x).*u)./H.(x)

    return uξ, uη, uσ
end
