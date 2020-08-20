using SparseArrays, PyPlot, LinearAlgebra

#close("all")
#plt.style.use("~/presentation_plots.mplstyle")

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

# parameters
L = 1e6
H0 = 1e3
Pr = 1
f = 1e-4
N = 1e-3
kap = 1e-1

# number of grid points (warning: use more than ~2^16 points at your own risk)
nx = 2^7
nz = 2^7

# domain (coordinates of cell centers)
dx = L/nx
x = dx/2:dx:L-dx/2
dz = H0/nz
z = -H0+dz/2:dz:-dz/2

println("Ekman layer thickness: ", sqrt(2*kap/f))
println("dz = ", dz)

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
        b[i, j] = N^2*1e1*sin(2*x[i]*pi/L)*exp(-(z[j] + H0)^2/2/(H0/4)^2)
        
        #= # gaussian in x and z (centered at center) =#
        #= b[i, j] = N^2*1e1*exp(-(z[j] + H0/2)^2/2/(H0/8)^2 - (x[i] - L/2)^2/2/(L/8)^2) =#
    end
end

"""
    A, rhs = assembleMatrix()

Setup linear system for problem.
"""
function assembleMatrix()
    nPts = nx*nz
    nVars = 4

    umap = reshape(1:nPts*nVars, nVars, nx, nz)    
    A = Tuple{Int64,Int64,Float64}[]  
    rhs = zeros(nPts*nVars)                      

    # Main loop, insert stencil in matrix for each node point
    for i = 1:nx, j = 1:nz

        # Boundary conditions: periodic in x
        iL = mod1(i-1, nx)
        iR = mod1(i+1, nx)

        # equation 1: -f*v + dx(p) - dz(nu*dz(u)) = 0
        row = umap[1, i, j] 
        # u
        if 1 < j < nz
          push!(A, (row, umap[1, i, j-1], -Pr*kap/dz^2))
          push!(A, (row, umap[1, i, j], 2Pr*kap/dz^2))
          push!(A, (row, umap[1, i, j+1], -Pr*kap/dz^2))
        elseif j == 1 # applying u = 0 BC
          push!(A, (row, umap[1, i, j], 3Pr*kap/dz^2))
          push!(A, (row, umap[1, i, j+1], -Pr*kap/dz^2))
        elseif j == nz # applying dz(u) = 0 BC
          push!(A, (row, umap[1, i, j-1], -Pr*kap/dz^2))
          push!(A, (row, umap[1, i, j], Pr*kap/dz^2))
        end
        # v
        push!(A, (row, umap[2, iL, j], -f/2))
        push!(A, (row, umap[2, i, j], -f/2))
        # p
        push!(A, (row, umap[4, iL, j], -1/dx))
        push!(A, (row, umap[4, i, j],  1/dx))
            
        # equation 2: f*u - dz(nu*dz(v)) = 0
        row = umap[2, i, j] 
        # v
        if 1 < j < nz
          push!(A, (row, umap[2, i, j-1], -Pr*kap/dz^2))
          push!(A, (row, umap[2, i, j], 2Pr*kap/dz^2))
          push!(A, (row, umap[2, i, j+1], -Pr*kap/dz^2))
        elseif j == 1 # applying v = 0 BC
          push!(A, (row, umap[2, i, j], 3Pr*kap/dz^2))
          push!(A, (row, umap[2, i, j+1], -Pr*kap/dz^2))
        elseif j == nz # applying dz(v) = 0 BC
          push!(A, (row, umap[2, i, j-1], -Pr*kap/dz^2))
          push!(A, (row, umap[2, i, j], Pr*kap/dz^2))
        end
        # u
        push!(A, (row, umap[1, iR, j], f/2))
        push!(A, (row, umap[1, i, j], f/2))

        # equation 3: dz(p) = b
        row = umap[3, i, j] 
        if j < nz
          # p
          push!(A, (row, umap[4, i, j], -1/dz))
          push!(A, (row, umap[4, i, j+1], 1/dz))
          # b
          rhs[row] = (b[i, j] + b[i, j+1])/2
        else # apply w = 0 at z = -H BC instead
          push!(A, (row, umap[3, i, 1], 1))
        end

        # equation 4: dx(u) + dz(w) = 0
        row = umap[4, i, j] 
        # u
        push!(A, (row, umap[1, i, j], -1/dx))
        push!(A, (row, umap[1, iR, j], 1/dx))
        # w
        push!(A, (row, umap[3, i, j], -1/dz))
        if j < nz # else: w = 0 at z = 0 from BC
          push!(A, (row, umap[3, i, j+1], 1/dz))
        end

    end

    # Create CSC sparse matrix from matrix elements
    A = sparse((x->x[1]).(A), (x->x[2]).(A), (x->x[3]).(A), nPts*nVars, nPts*nVars)

    println("rank(A)    = ", rank(A))
    println("nPts*nVars = ", nPts*nVars)

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
    sol = reshape(A \ rhs, 4, nx, nz)

    return sol
end

"""
    make_plots(sol)

Plot the results of solver.
"""
function make_plots(sol)
    # unpack variables
    u = sol[1, :, :]
    v = sol[2, :, :]
    w = sol[3, :, :]
    p = sol[4, :, :]
    rho = N^2*repeat(z', nx, 1) + b 

    # coordinates of cell corners
    xu = -dx/2:dx:L-dx/2
    zu = -H0:dz:0

    xv = 0:dx:L
    zv = -H0:dz:0

    xw = -dx/2:dx:L-dx/2
    zw = -H0-dz/2:dz:-dz/2

    xp = 0:dx:L
    zp = -H0:dz:0

    # plot each
    figure()
    umax = maximum(abs.(u))
    umin = -umax
    pcolormesh(xu, zu, u', cmap="RdBu_r", vmin=umin, vmax=umax)
    colorbar(label=L"$u$ (m s$^{-1}$)")
    contour(x, z, rho', 10, colors="k", alpha=0.5)
    title(L"$u$")
    xlabel(L"$x$ (m)")
    ylabel(L"$z$ (m)")

    figure()
    vmax = maximum(abs.(v))
    vmin = -vmax
    pcolormesh(xv, zv, v', cmap="RdBu_r", vmin=vmin, vmax=vmax)
    colorbar(label=L"$v$ (m s$^{-1}$)")
    contour(x, z, rho', 10, colors="k", alpha=0.5)
    title(L"$v$")
    xlabel(L"$x$ (m)")
    ylabel(L"$z$ (m)")

    figure()
    vmax = maximum(abs.(w))
    vmin = -vmax
    pcolormesh(xw, zw, w', cmap="RdBu_r", vmin=vmin, vmax=vmax)
    colorbar(label=L"$w$ (m s$^{-1}$)")
    contour(x, z, rho', 10, colors="k", alpha=0.5)
    title(L"$w$")
    xlabel(L"$x$ (m)")
    ylabel(L"$z$ (m)")

    figure()
    pcolormesh(xp, zp, p')
    #= plot(x[2], z[2], "ro") # to show where p = 0 is set =#
    colorbar(label=L"$p$ (m$^2$ s$^{-2}$)")
    contour(x, z, rho', 10, colors="k", alpha=0.5)
    title(L"$p$")
    xlabel(L"$x$ (m)")
    ylabel(L"$z$ (m)")
end


"""
    ux, wz = check_continuity(sol)

See how true continuity equation ux + wz = 0 is and plot results.
"""
function check_continuity(sol)
    # unpack variables
    u = sol[1, :, :]
    v = sol[2, :, :]
    w = sol[3, :, :]
    p = sol[4, :, :]
    rho = N^2*repeat(z', nx, 1) + b 
    
    # compute dx(u)
    ux = zeros(nx, nz)
    for i=1:nx-1
        ux[i, :] = (u[i+1, :] - u[i, :])/dx
    end
    ux[nx, :] = (u[1, :] - u[nx, :])/dx

    # compute dz(w)
    wz = zeros(nx, nz)
    for j=1:nz-1
      wz[:, j] = (w[:, j+1] - w[:, j])/dz
    end
    wz[:, nz] = -w[:, nz]/dz

    # plot and print max dx(u) + dz(w)
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

# main
sol = solve()
make_plots(sol)
ux, wz = check_continuity(sol)
