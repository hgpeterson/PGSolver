"""
    c = mkfdstencil(x, xbar, k)

Compute the coefficients `c` in a finite difference approximation of a function
defined at the grid points `x`, evaluated at `xbar`, of order `k`.
"""
function mkfdstencil(x, xbar, k)
	n = length(x)
	A = @. (x[:]' - xbar) ^ (0:n-1) / factorial(0:n-1)
	b = zeros(n)
	b[k+1] = 1.0
	c = A \ b
end

"""
    integral = trapz(f, x)

Integrate array `f` over domain `x` using trapezoidal rule.
"""
function trapz(f, x)
    return 0.5*sum((f[1:end-1] .+ f[2:end]).*(x[2:end] .- x[1:end-1]))
end

"""
    integral = cumtrapz(f, x)

Cumulatively integrate array `f` over domain `x` using trapezoidal rule.
"""
function cumtrapz(f, x)
    y = zeros(size(f, 1))
    for i=2:size(f, 1)
        y[i] = trapz(f[1:i], x[1:i])
    end
    return y
end

"""
    fz = differentiate(f, z)

Compute second order first derivative of `f` on grid `z`.
"""
function differentiate(f, z)
    # allocate derivative array, fz
    nz = size(z, 1)
    fz = zeros(size(f))

    # 2nd order centered difference
    for j=2:nz-1
        fd_z = mkfdstencil(z[j-1:j+1], z[j], 1)
        fz[j] = sum(fd_z.*f[j-1:j+1])
    end

    # 2nd order off-center difference on top and bottom boundary
    fd_bot = mkfdstencil(z[1:3], z[1], 1)
    fd_top = mkfdstencil(z[nz-2:nz], z[nz], 1)
    fz[1]  = sum(fd_bot.*f[1:3])
    fz[nz] = sum(fd_top.*f[nz-2:nz])

    return fz
end

"""
    fr = rDerivativeTF(field)

Compute dr(`field`) in terrian-following coordinates.
Note: dr() = dρ() - dr(H)*σ*dσ()/H
"""
function rDerivativeTF(field)
    # allocate
    fr = zeros(nρ, nσ)

    # dρ(field)
    for j=1:nσ
        fr[:, j] += differentiate(field[:, j], ρ)
    end

    # -dr(H)*σ*dσ(field)/H
    for i=1:nρ
        fr[i, :] -= Hr(ρ[i])*σ.*differentiate(field[i, :], σ)/H(ρ[i])
    end

    return fr
end

"""
    fz = zDerivativeTF(field)

Compute dz(`field`) in terrian-following coordinates.
Note: dz() = dσ()/H
"""
function zDerivativeTF(field)
    # allocate
    fz = zeros(nρ, nσ)

    # dσ(field)/H
    for i=1:nρ
        fz[i, :] += differentiate(field[i, :], σ)/H(ρ[i])
    end

    return fz
end
