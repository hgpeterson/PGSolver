using SparseArrays, PyPlot, LinearAlgebra

close("all")
plt.style.use("~/presentation_plots.mplstyle")

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

# parameters
L = 1e6
H0 = 1e3
Pr = 1
f = 1e-4
N = 1e-3

# number of grid points
nx = 2^7 + 1 
nz = 2^7

# domain
dx = L/(nx - 1)
x = 0:dx:L
dz = H0/(nz - 1)
z = -H0:dz:0
#= z = @. -H0*(cos(pi*(0:nz-1)/(nz-1)) + 1)/2 # chebyshev =# 

# diffusivity
kap = 1e0 # (constant for now for easier stencil)
println("Ekman layer thickness: ", sqrt(2*kap/f))
println("z[2] - z[1] = ", z[2] - z[1])

# buoyancy perturbation
b = zeros(nx, nz)
for i=1:nx
    for j=1:nz
        # constant
        #= b[i, j] = 1 =#

        # gaussian in z (centered at bottom)
        #= b[i, j] = N^2*1e2*exp(-(z[j] + H0)^2/2/(H0/4)^2) =#
        
        # sine in x
        #= b[i, j] = N^2*1e1*sin(2*x[i]*pi/(L + dx)) =#

        # sine in x, gaussian in z (centered at bottom)
        b[i, j] = N^2*1e1*sin(2*x[i]*pi/(L + dx))*exp(-(z[j] + H0)^2/2/(H0/4)^2)
        
        # gaussian in x and z (centered at center)
        #= b[i, j] = N^2*1e1*exp(-(z[j] + H0/2)^2/2/(H0/8)^2 - (x[i] - L/2)^2/2/(L/8)^2) =#
    end
end

"""
    A, rhs = assembleMatrix()

Setup linear system for problem.
"""
function assembleMatrix()
    nPts = nx*nz
    iU = nPts + 1 # add equation for vertically integrated zonal flow

    umap = reshape(1:nPts, nx, nz)    
    A = Tuple{Int64,Int64,Float64}[]  
    rhs = zeros(nPts + 1)                      

    # for finite difference on the top and bottom boundary
    fd_bot = mkfdstencil(z[1:3], z[1], 1)
    fd_top = mkfdstencil(z[nz-2:nz], z[nz], 1)
    fd_top_zz = mkfdstencil(z[nz-3:nz], z[nz], 2)

    # Main loop, insert stencil in matrix for each node point
    for i=1:nx
        # Boundary conditions: periodic in x
        iL = mod1(i-1, nx)
        iR = mod1(i+1, nx)

        # Lower boundary conditions 
        # b.c. 1: dz(chi) = 0
        push!(A, (umap[i, 1], umap[i, 1], fd_bot[1]))
        push!(A, (umap[i, 1], umap[i, 2], fd_bot[2]))
        push!(A, (umap[i, 1], umap[i, 3], fd_bot[3]))
        # b.c. 2: chi = 0 
        push!(A, (umap[i, 2], umap[i, 1], 1.0))

        # Upper boundary conditions
        # b.c. 1: dzz(chi) = 0 (or -stress at top)
        push!(A, (umap[i, nz], umap[i, nz-3], fd_top_zz[1]))
        push!(A, (umap[i, nz], umap[i, nz-2], fd_top_zz[2]))
        push!(A, (umap[i, nz], umap[i, nz-1], fd_top_zz[3]))
        push!(A, (umap[i, nz], umap[i, nz],   fd_top_zz[4]))
        #= rhs[umap[1, i, nz]] = -0.1 # some wind stress at the top =#
        # b.c. 2: chi = bottom stress 
        push!(A, (umap[i, nz-1], umap[i, nz], 1.0))
        push!(A, (umap[i, nz-1], iU, -1.0))

        # Interior nodes
        for j=3:nz-2
            row = umap[i, j] 

            # dz stencil
            fd_z = mkfdstencil([z[j-1] z[j] z[j+1]], z[j], 1)

            # dzz stencil
            fd_zz = mkfdstencil([z[j-1] z[j] z[j+1]], z[j], 2)

            # dzzzz stencil
            fd_zzzz = mkfdstencil(z[j-2:j+2], z[j], 4)
            
            # eqtn: dzz(nu*dzz(chi)) + f^2*chi/nu + f*U = dx(b)
            # term 1
            push!(A, (row, umap[i, j-2], Pr*kap*fd_zzzz[1]))
            push!(A, (row, umap[i, j-1], Pr*kap*fd_zzzz[2]))
            push!(A, (row, umap[i, j],   Pr*kap*fd_zzzz[3]))
            push!(A, (row, umap[i, j+1], Pr*kap*fd_zzzz[4]))
            push!(A, (row, umap[i, j+2], Pr*kap*fd_zzzz[5]))
            # term 2
            push!(A, (row, umap[i, j], f^2/(Pr*kap)))
            # term 3
            push!(A, (row, iU, f))
            # term 4 (rhs)
            rhs[row] = (b[iR, j] - b[iL, j])/(2*dx)
        end
    end

    # zonal mean / integral equation for bottom stress
    #   < int_-H^0 (f^2*chi/nu + U) dz + dz(nu*dzz(chi)) > = 0
    row = iU
    #= push!(A, (row, iU, 1.0)) =#
    # term 1
    for i=1:nx
        for j=1:nz
            push!(A, (row, umap[i, j], f^2/(Pr*kap)*dx/L*dz))
            push!(A, (row, iU,   f*dx/L*dz))
        end
    end
    # term 2
    fd_zzz_top = mkfdstencil(z[nz-4:nz], z[nz], 3)
    for i=1:nx
        push!(A, (row, umap[i, nz-4], Pr*kap*fd_zzz_top[1]*dx/L))
        push!(A, (row, umap[i, nz-3], Pr*kap*fd_zzz_top[2]*dx/L))
        push!(A, (row, umap[i, nz-2], Pr*kap*fd_zzz_top[3]*dx/L))
        push!(A, (row, umap[i, nz-1], Pr*kap*fd_zzz_top[4]*dx/L))
        push!(A, (row, umap[i, nz],   Pr*kap*fd_zzz_top[5]*dx/L))
    end

    # Create CSC sparse matrix from matrix elements
    A = sparse((x->x[1]).(A), (x->x[2]).(A), (x->x[3]).(A), nPts + 1, nPts + 1)

    #= println("rank(A)  = ", rank(A)) =#
    #= println("nPts + 1 = ", nPts + 1) =#

    return A, rhs
end

"""
    sol = solve()

Setup and solve linear system.
"""
function solve()
    # get matrix and vector
    A, rhs = assembleMatrix()

    # Solve + reshape for plotting
    sol = A \ rhs
    println("U = ", sol[end])
    sol = reshape(sol[1:end-1], nx, nz)

    return sol
end

"""
    u, w = make_plots(sol)

Plot the results of solver.
"""
function make_plots(sol)
    # gather solution
    chi = sol[:, :]
    rho = N^2*repeat(z', nx, 1) + b 

    # compute u = dz(chi)
    u = zeros(nx, nz)
    for j=2:nz-1
        fd_z = mkfdstencil(z[j-1:j+1], z[j], 1)
        u[:, j] = (fd_z[1]*chi[:, j-1] + fd_z[2]*chi[:, j] + fd_z[3]*chi[:, j+1])
    end
    fd_z = mkfdstencil(z[1:3], z[1], 1)
    u[:, 1] = (fd_z[1]*chi[:, 1] + fd_z[2]*chi[:, 2] + fd_z[3]*chi[:, 3])
    fd_z = mkfdstencil(z[nz-2:nz], z[nz], 1)
    u[:, nz] = (fd_z[1]*chi[:, nz-2] + fd_z[2]*chi[:, nz-1] + fd_z[3]*chi[:, nz])

    # compute v = int(f*chi/nu)
    v = zeros(nx, nz)
    for i=1:nx
        v[i, :] = cumtrapz(f*chi[i, :]/(Pr*kap), z)
    end

    # compute w = -dx(chi)
    w = zeros(nx, nz)
    for i=2:nx-1
        w[i, :] = -(chi[i+1, :] - chi[i-1, :])/(2*dx)
    end
    w[1, :] = -(chi[2, :] - chi[nx, :])/(2*dx)
    w[nx, :] = -(chi[1, :] - chi[nx-1, :])/(2*dx)

    # plot
    figure()
    vmax = maximum(abs.(chi))
    vmin = -vmax
    pcolormesh(x, z, chi', cmap="RdBu_r", vmin=vmin, vmax=vmax)
    colorbar(label=L"$\chi$ (m$^2$ s$^{-1}$)")
    contour(x, z, rho', 10, colors="k", alpha=0.5)
    title(L"$\chi$")
    xlabel(L"$x$ (m)")
    ylabel(L"$z$ (m)")

    figure()
    vmax = maximum(abs.(u))
    vmin = -vmax
    pcolormesh(x, z, u', cmap="RdBu_r", vmin=vmin, vmax=vmax)
    colorbar(label=L"$u$ (m s$^{-1}$)")
    contour(x, z, rho', 10, colors="k", alpha=0.5)
    title(L"$u$")
    xlabel(L"$x$ (m)")
    ylabel(L"$z$ (m)")

    figure()
    vmax = maximum(abs.(v))
    vmin = -vmax
    pcolormesh(x, z, v', cmap="RdBu_r", vmin=vmin, vmax=vmax)
    colorbar(label=L"$v$ (m s$^{-1}$)")
    contour(x, z, rho', 10, colors="k", alpha=0.5)
    title(L"$v$")
    xlabel(L"$x$ (m)")
    ylabel(L"$z$ (m)")

    figure()
    vmax = maximum(abs.(w))
    vmin = -vmax
    pcolormesh(x, z, w', cmap="RdBu_r", vmin=vmin, vmax=vmax)
    colorbar(label=L"$w$ (m s$^{-1}$)")
    contour(x, z, rho', 10, colors="k", alpha=0.5)
    title(L"$w$")
    xlabel(L"$x$ (m)")
    ylabel(L"$z$ (m)")

    return u, w
end

"""
    ux, wz = check_continuity(sol, u, w)

See how true continuity equation ux + wz = 0 is and plot results.
"""
function check_continuity(sol, u, w)
    # gather solution
    chi = sol[:, :]
    rho = N^2*repeat(z', nx, 1) + b 

    # compute dx(u)
    ux = zeros(nx, nz)
    for i=2:nx-1
        ux[i, :] = (u[i+1, :] - u[i-1, :])/(2*dx)
    end
    ux[1, :] = (u[2, :] - u[nx, :])/(2*dx)
    ux[nx, :] = (u[1, :] - u[nx-1, :])/(2*dx)

    # compute dz(w)
    wz = zeros(nx, nz)
    for j=2:nz-1
        fd_z = mkfdstencil(z[j-1:j+1], z[j], 1)
        wz[:, j] = fd_z[1]*w[:, j-1] + fd_z[2]*w[:, j] + fd_z[3]*w[:, j+1]
    end
    fd_z = mkfdstencil(z[1:3], z[1], 1)
    wz[:, 1] = fd_z[1]*w[:, 1] + fd_z[2]*w[:, 2] + fd_z[3]*w[:, 3]
    fd_z = mkfdstencil(z[nz-2:nz], z[nz], 1)
    wz[:, nz] = fd_z[1]*w[:, nz-2] + fd_z[2]*w[:, nz-1] + fd_z[3]*w[:, nz]

    # plot and print maximum of dx(u) + dz(w)
    figure()
    vmax = maximum(abs.(ux))
    vmin = -vmax
    pcolormesh(x, z, ux', cmap="RdBu_r", vmin=vmin, vmax=vmax)
    colorbar(label=L"$u_x$ (m s$^{-1}$)")
    contour(x, z, rho', 10, colors="k", alpha=0.5)
    title(L"$u_x$")
    xlabel(L"$x$ (m)")
    ylabel(L"$z$ (m)")

    figure()
    vmax = maximum(abs.(wz))
    vmin = -vmax
    pcolormesh(x, z, wz', cmap="RdBu_r", vmin=vmin, vmax=vmax)
    colorbar(label=L"$w_z$ (m s$^{-1}$)")
    contour(x, z, rho', 10, colors="k", alpha=0.5)
    title(L"$w_z$")
    xlabel(L"$x$ (m)")
    ylabel(L"$z$ (m)")

    figure()
    vmax = maximum(abs.(ux .+ wz))
    ijmax = argmax(abs.(ux .+ wz))
    println("Max ux + wz = ", vmax, " at ", ijmax)
    vmin = -vmax
    pcolormesh(x, z, (ux .+ wz)', cmap="RdBu_r", vmin=vmin, vmax=vmax)
    colorbar(label=L"$u_x + w_z$ (m s$^{-1}$)")
    contour(x, z, rho', 10, colors="k", alpha=0.5)
    title(L"$u_x + w_z$")
    xlabel(L"$x$ (m)")
    ylabel(L"$z$ (m)")

    return ux, wz
end

sol = solve()
u, w = make_plots(sol)
#= check_continuity(sol, u, w) =#
