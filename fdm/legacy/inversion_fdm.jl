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

# parameters
L = 1e6
H0 = 1e3
Pr = 1
f = 1e-4
N = 1e-3

# number of grid points (warning: use more than ~2^16 points at your own risk)
nx = 2^7 + 1 
nz = 2^7

# domain
dx = L/(nx - 1)
x = 0:dx:L
#= dz = H0/(nz - 1) =#
#= z = -H0:dz:0 =#
z = @. -H0*(cos(pi*(0:nz-1)/(nz-1)) + 1)/2 # chebyshev 

# diffusivity
kap = 1e0*ones(nz)
#= kap = @. 1e-1 + 1e1*exp(-(z + H0)/2e2) =#
println("Ekman layer thickness: ", sqrt(2*kap[1]/f))
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

    # for finite difference on the top and bottom boundary
    fd_bot = mkfdstencil(z[1:3], z[1], 1)
    fd_top = mkfdstencil(z[nz-2:nz], z[nz], 1)

    # Main loop, insert stencil in matrix for each node point
    for i = 1:nx
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

        # Lower boundary conditions (j = 1)
        # b.c. 1: u = 0
        push!(A, (umap[1, i, 1], umap[1, i, 1], 1.0))
        # b.c. 2: v = 0
        push!(A, (umap[2, i, 1], umap[2, i, 1], 1.0))
        # b.c. 3: w = 0
        push!(A, (umap[3, i, 1], umap[3, i, 1], 1.0))
        # b.c. 4: dz(p) = b 
        # p
        push!(A, (umap[4, i, 1], umap[4, i, 1], fd_bot[1]))
        push!(A, (umap[4, i, 1], umap[4, i, 2], fd_bot[2]))
        push!(A, (umap[4, i, 1], umap[4, i, 3], fd_bot[3]))
        # b
        rhs[umap[4, i, 1]] = b[i, 1]

        # Upper boundary conditions (j = nz)
        # b.c. 1: dz(u) = 0
        push!(A, (umap[1, i, nz], umap[1, i, nz-2], fd_top[1]))
        push!(A, (umap[1, i, nz], umap[1, i, nz-1], fd_top[2]))
        push!(A, (umap[1, i, nz], umap[1, i, nz],   fd_top[3]))
        #= rhs[umap[1, i, nz]] = 0.1 # if you want to add wind stress =#
        # b.c. 2: dz(v) = 0
        push!(A, (umap[2, i, nz], umap[2, i, nz-2], fd_top[1]))
        push!(A, (umap[2, i, nz], umap[2, i, nz-1], fd_top[2]))
        push!(A, (umap[2, i, nz], umap[2, i, nz],   fd_top[3]))
        #= rhs[umap[2, i, nz]] = 0.1 # if you want to add wind stress =#
        # b.c. 3: w = 0
        push!(A, (umap[3, i, nz], umap[3, i, nz], 1.0))
        # b.c. 4: dz(p) = b
        # p
        push!(A, (umap[4, i, nz], umap[4, i, nz-2], fd_top[1]))
        push!(A, (umap[4, i, nz], umap[4, i, nz-1], fd_top[2]))
        push!(A, (umap[4, i, nz], umap[4, i, nz],   fd_top[3]))
        # b
        rhs[umap[4, i, nz]] = b[i, nz]

        ###
        # begin 
        #
        # testing where to set p = 0:
        #   replace relevant lines above with desired code below
        ###
        
        # note: replace a line like
        #
        #     push!(A, (umap[?, i, j], umap[4, i, j], 1.0)) 
        #
        # with
        #
        #     for k=1:nx
        #         for l=1:nz
        #             push!(A, (umap[?, i, j], umap[4, k, l], 1.0))
        #         end
        #     end
        #
        # if you want to set sum(p) = 0 instead of p = 0 at a point.
        
        # chosen point
        #= iSwap = 20 =#

        #= # instead of w = 0 at bottom (result: singular) =#
        #= if i == iSwap =#
        #=     push!(A, (umap[3, i, 1], umap[4, i, 1], 1.0)) =#
        #= else =#
        #=     push!(A, (umap[3, i, 1], umap[3, i, 1], 1.0)) =#
        #= end =#

        # instead of w = 0 at top (result: singular)
        #= if i == iSwap =#
        #=     push!(A, (umap[3, i, nz], umap[4, i, nz], 1.0)) =#
        #= else =#
        #=     push!(A, (umap[3, i, nz], umap[3, i, nz], 1.0)) =#
        #= end =#
        
        #= # instead of dz(p) = b at top (result: leaves a messed up column) =#
        #= if i == iSwap =#
        #=     push!(A, (umap[4, i, nz], umap[4, i, nz], 1.0)) =#
        #= else =#
        #=     # p =#
        #=     push!(A, (umap[4, i, nz], umap[4, i, nz-2], fd_top[1])) =#
        #=     push!(A, (umap[4, i, nz], umap[4, i, nz-1], fd_top[2])) =#
        #=     push!(A, (umap[4, i, nz], umap[4, i, nz],   fd_top[3])) =#
        #=     # b =#
        #=     rhs[umap[4, i, nz]] = b[i, nz] =#
        #= end =#

        ###
        # end
        ###
        

        for j = 2:nz-1
            # Interior nodes

            # dz stencil
            fd_z = mkfdstencil([z[j-1] z[j] z[j+1]], z[j], 1)

            # dzz stencil
            fd_zz = mkfdstencil([z[j-1] z[j] z[j+1]], z[j], 2)
            
            # dz(kap) at z = z[j]
            dzkap = fd_z[1]*kap[j-1] + fd_z[2]*kap[j] + fd_z[3]*kap[j+1]

            # equation 1: -f*v + dx(p) - dz(nu*dz(u)) = 0
            row = umap[1, i, j] 
            # u
            push!(A, (row, umap[1, i, j-1], -Pr*(kap[j]*fd_zz[1] + fd_z[1]*dzkap)))
            push!(A, (row, umap[1, i, j],   -Pr*(kap[j]*fd_zz[2] + fd_z[2]*dzkap)))
            push!(A, (row, umap[1, i, j+1], -Pr*(kap[j]*fd_zz[3] + fd_z[3]*dzkap)))
            # v
            push!(A, (row, umap[2, i, j], -f))
            # p
            push!(A, (row, umap[4, iL, j], -1.0/(2*dx)))
            push!(A, (row, umap[4, iR, j],  1.0/(2*dx)))
            
            # equation 2: f*u - dz(nu*dz(v)) = 0
            row = umap[2, i, j] 
            # u
            push!(A, (row, umap[1, i, j], f))
            # v
            push!(A, (row, umap[2, i, j-1], -Pr*(kap[j]*fd_zz[1] + fd_z[1]*dzkap)))
            push!(A, (row, umap[2, i, j],   -Pr*(kap[j]*fd_zz[2] + fd_z[2]*dzkap)))
            push!(A, (row, umap[2, i, j+1], -Pr*(kap[j]*fd_zz[3] + fd_z[3]*dzkap)))

            # equation 3: dx(u) + dz(w) = 0
            row = umap[3, i, j] 
            # u
            push!(A, (row, umap[1, iL, j], -1.0/(2*dx)))
            push!(A, (row, umap[1, iR, j],  1.0/(2*dx)))
            # w
            push!(A, (row, umap[3, i, j-1], fd_z[1]))
            push!(A, (row, umap[3, i, j],   fd_z[2]))
            push!(A, (row, umap[3, i, j+1], fd_z[3]))
            
            ###
            # if you wish to replace continuity at a point with p = 0, use 
            # this instead of the above block
            #   result: singular
            ###
            #= row = umap[3, i, j] =# 
            #= if i == 2 && j == 2 =#
            #=     push!(A, (row, umap[4, i, j], 1.0)) =#
            #= else =#
            #=     # u =#
            #=     push!(A, (row, umap[1, iL, j], -1.0/(2*dx))) =#
            #=     push!(A, (row, umap[1, iR, j],  1.0/(2*dx))) =#
            #=     # w =#
            #=     push!(A, (row, umap[3, i, j-1], fd_z[1])) =#
            #=     push!(A, (row, umap[3, i, j],   fd_z[2])) =#
            #=     push!(A, (row, umap[3, i, j+1], fd_z[3])) =#
            #= end =#
            ###

            # equation 4: dz(p) = b
            row = umap[4, i, j] 
            # p
            push!(A, (row, umap[4, i, j-1], fd_z[1]))
            push!(A, (row, umap[4, i, j],   fd_z[2]))
            push!(A, (row, umap[4, i, j+1], fd_z[3]))
            # b
            rhs[row] = b[i, j]
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

    # plot each
    figure()
    umax = maximum(abs.(u))
    umin = -umax
    pcolormesh(x, z, u', cmap="bwr", vmin=umin, vmax=umax)
    colorbar(label=L"$u$ (m s$^{-1}$)")
    contour(x, z, rho', 10, colors="k", alpha=0.5)
    title(L"$u$")
    xlabel(L"$x$ (m)")
    ylabel(L"$z$ (m)")

    figure()
    vmax = maximum(abs.(v))
    vmin = -vmax
    pcolormesh(x, z, v', cmap="bwr", vmin=vmin, vmax=vmax)
    colorbar(label=L"$v$ (m s$^{-1}$)")
    contour(x, z, rho', 10, colors="k", alpha=0.5)
    title(L"$v$")
    xlabel(L"$x$ (m)")
    ylabel(L"$z$ (m)")

    figure()
    vmax = maximum(abs.(w))
    vmin = -vmax
    pcolormesh(x, z, w', cmap="bwr", vmin=vmin, vmax=vmax)
    colorbar(label=L"$w$ (m s$^{-1}$)")
    contour(x, z, rho', 10, colors="k", alpha=0.5)
    title(L"$w$")
    xlabel(L"$x$ (m)")
    ylabel(L"$z$ (m)")

    figure()
    pcolormesh(x, z, p')
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

    # plot and print max dx(u) + dz(w)
    figure()
    vmax = maximum(abs.(ux))
    vmin = -vmax
    pcolormesh(x, z, ux', cmap="bwr", vmin=vmin, vmax=vmax)
    colorbar(label=L"$u_x$ (m s$^{-1}$)")
    contour(x, z, rho', 10, colors="k", alpha=0.5)
    title(L"$u_x$")
    xlabel(L"$x$ (m)")
    ylabel(L"$z$ (m)")

    figure()
    vmax = maximum(abs.(wz))
    vmin = -vmax
    pcolormesh(x, z, wz', cmap="bwr", vmin=vmin, vmax=vmax)
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
    pcolormesh(x, z, (ux .+ wz)', cmap="bwr", vmin=vmin, vmax=vmax)
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
