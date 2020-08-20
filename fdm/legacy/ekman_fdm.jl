# Solve boussinesq equations on a slope on a square (y, z)-grid using the 
# finite difference method.

using SparseArrays, PyPlot, LinearAlgebra

#= close("all") =#
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

# parameters
L = 1e6
H0 = 1e3
Pr = 1
f = 1e-4

# number of grid points
nx = 2^7
nz = 2^7

# domain
x = 0:L/(nx-1):L
#= z = -H0:H0/(nz-1):0 =#
z = @. -H0*(cos(pi*(0:nz-1)/(nz-1)) + 1)/2 # chebyshev nodes

# coeffs
kap = 1e0*ones(nz)

function assembleMatrix()
    N = 2*nx*nz

    umap = reshape(1:N, 2, nx, nz)    
    A = Tuple{Int64,Int64,Float64}[]  
    b = zeros(N)                      

    # for finite difference on the top and bottom boundary
    fd_bot = mkfdstencil(z[1:3], z[1], 1)
    fd_top = mkfdstencil(z[nz-2:nz], z[nz], 1)

    # Main loop, insert stencil in matrix for each node point
    for i = 1:nx
        for j = 1:nz
            if j == 1 
                # Lower boundary conditions
                # b.c. 1: u = u_I
                push!(A, (umap[1, i, j], umap[1, i, j], 1.0))
                b[umap[1, i, j]] = 1e-2
                # b.c. 2: v = v_I
                push!(A, (umap[2, i, j], umap[2, i, j], 1.0))
                b[umap[2, i, j]] = 1e-2
            elseif j == nz
                # Upper boundary conditions
                
                # b.c. 1: u_z = 0
                push!(A, (umap[1, i, j], umap[1, i, j-2], fd_top[1]))
                push!(A, (umap[1, i, j], umap[1, i, j-1], fd_top[2]))
                push!(A, (umap[1, i, j], umap[1, i, j],   fd_top[3]))
                # b.c. 2: v_z = 0
                push!(A, (umap[2, i, j], umap[2, i, j-2], fd_top[1]))
                push!(A, (umap[2, i, j], umap[2, i, j-1], fd_top[2]))
                push!(A, (umap[2, i, j], umap[2, i, j],   fd_top[3]))
            else
                # Interior nodes
               
                # Boundary conditions: periodic in x
                if i == 1 
                    iL = nx
                    iR = 2
                elseif i == nx
                    iL = nx - 1
                    iR = 1
                else
                    iL = i - 1
                    iR = i + 1
                end

                # dz stencil
                fd_z = mkfdstencil([z[j-1] z[j] z[j+1]], z[j], 1)

                # dzz stencil
                fd_zz = mkfdstencil([z[j-1] z[j] z[j+1]], z[j], 2)
                
                # dz(kap) at z = z[j]
                dzkap = fd_z[1]*kap[j-1] + fd_z[2]*kap[j] + fd_z[3]*kap[j+1]

                # equation 1: -f*v - dz(nu*dz(u)) = 0
                row = umap[1, i, j] 
                # u
                push!(A, (row, umap[1, i, j-1], -Pr*(kap[j]*fd_zz[1] + fd_z[1]*dzkap)))
                push!(A, (row, umap[1, i, j],   -Pr*(kap[j]*fd_zz[2] + fd_z[2]*dzkap)))
                push!(A, (row, umap[1, i, j+1], -Pr*(kap[j]*fd_zz[3] + fd_z[3]*dzkap)))
                # v
                push!(A, (row, umap[2, i, j], -f))
                
                # equation 2: f*u - dz(nu*dz(v)) = 0
                row = umap[2, i, j] 
                # u
                push!(A, (row, umap[1, i, j], f))
                # v
                push!(A, (row, umap[2, i, j-1], -Pr*(kap[j]*fd_zz[1] + fd_z[1]*dzkap)))
                push!(A, (row, umap[2, i, j],   -Pr*(kap[j]*fd_zz[2] + fd_z[2]*dzkap)))
                push!(A, (row, umap[2, i, j+1], -Pr*(kap[j]*fd_zz[3] + fd_z[3]*dzkap)))
            end
        end
    end

    # Create CSC sparse matrix from matrix elements
    A = sparse((x->x[1]).(A), (x->x[2]).(A), (x->x[3]).(A), N, N)

    return A, b
end

function solve()
    A, b = assembleMatrix()

    # Solve + reshape for plotting
    sol = reshape(A \ b, 2, nx, nz)

    return sol
end

function make_plots(sol)
    u = sol[1, :, :]
    v = sol[2, :, :]

    figure()
    umax = maximum(abs.(u))
    umin = -umax
    pcolormesh(x, z, u', cmap="bwr", vmin=umin, vmax=umax)
    colorbar(label=L"$u$ (m s$^{-1}$)")
    title(L"$u$")
    xlabel(L"$y$ (m)")
    ylabel(L"$z$ (m)")

    figure()
    vmax = maximum(abs.(v))
    vmin = -vmax
    pcolormesh(x, z, v', cmap="bwr", vmin=vmin, vmax=vmax)
    colorbar(label=L"$v$ (m s$^{-1}$)")
    title(L"$v$")
    xlabel(L"$y$ (m)")
    ylabel(L"$z$ (m)")
end

sol = solve()
make_plots(sol)
