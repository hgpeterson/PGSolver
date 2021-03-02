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
    fz = differentiate_pointwise(f, z, z0, n)

Compute `n`th order derivative of `f` at `z0` given grid `z`.
"""
function differentiate_pointwise(f, z, z0, n)
    fd_z = mkfdstencil(z, z0, n)
    return sum(fd_z.*f)
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
        fz[j] = differentiate_pointwise(f[j-1:j+1], z[j-1:j+1], z[j], 1)
    end

    # 2nd order off-center difference on top and bottom boundary
    fz[1] = differentiate_pointwise(f[1:3], z[1:3], z[1], 1)
    fz[nz] = differentiate_pointwise(f[nz-2:nz], z[nz-2:nz], z[nz], 1)

    return fz
end

"""
    dy = RK4(t, Δt, y, f)

Compute `dy`, the change in `y` given a timestep of `Δt` and that `dt(y, t) = f(y, t)`.
"""
function RK4(t, Δt, y, f)
	f1 = f(y, t)
    f2 = f(y + Δt*f1/2, t + Δt/2)
    f3 = f(y + Δt*f2/2, t + Δt/2)
    f4 = f(y + Δt*f3, t + Δt)
	dy = Δt/6*(f1 + 2*f2 + 2*f3 + f4)
    return dy
end
