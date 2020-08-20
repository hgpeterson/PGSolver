using SparseArrays, PyPlot, LinearAlgebra
include("setParams.jl")
include("mesh.jl")

close("all")
plt.style.use("~/presentation_plots.mplstyle")

"""
compute area of triangle defined by pts
"""
function tri_area(pts)
	p1 = pts[1, :]; p2 = pts[2, :]; p3 = pts[3, :]
	area = 1/2 * abs(p1[1]*(p2[2]-p3[2]) + p2[1]*(p3[2]-p1[2]) + p3[1]*(p1[2]-p2[2]))
    return area
end

"""
second order gaussian quadrature on a line
"""
function gaussian_quad2(func, a, b)
    # integration points 
    integration_pts = (b - a)/2 .*[-1/sqrt(3) 1/sqrt(3)] .+ (a + b)/2
    
    # weights are constant
    weights = [1 1]

    # return sum
    integral = 0
    for i=1:2
        integral += weights[i]*func(integration_pts[i])
    end

    return (b - a)/2*integral
end

"""
second order gaussian quadrature on a triangle
"""
function gaussian_quad2_tri(func, pts)
    # area of triangle
    area = tri_area(pts)

    # integration points (rows: point number, columns: x, z)
	integration_pts = 1/6*sum(pts, dims=1) .+ 1/2*I*pts
    
    # weights are constant
    weights = [1/3 1/3 1/3]

    # return sum
    integral = 0
    for i=1:3
        integral += weights[i]*func(integration_pts[i, :]')
    end

    return area*integral
end

"""
evaluate local basis function defined by cs = [c1 c2 c3] at pt = [x z]
"""
function local_basis_func(cs, pt)
    value = cs*[1 pt]'
    return value[1]
end

"""
	sol = fem_inversion(p, t, e, bx)
Inverts buoyancy field `b` for flow field on domain described by unstructured
triangular mesh `p, t`. 
"""
function fem_inversion(p, t, e, bx)
	nPts = size(p, 1)
	nTri = size(t, 1)

    nVars = 2

    umap = reshape(1:(nPts*nVars), nVars, nPts) 

	# Create global linear system using stamping method
    A = Tuple{Int64,Int64,Float64}[]  
    rhs = zeros(nPts*nVars)
	for k = 1:nTri
		# Calculate Elementary matrices and loads

		## get coeffs for linear basis func c1 + c2 x + c3 z
		V = zeros(3, 3)
		for row = 1:3
			V[row, :] = [1 p[t[k, row], 1] p[t[k, row], 2]]
		end
		C = inv(V)

		## calculate entries of matrix
        M_elem = zeros(3, 3)
        for i=1:3
            for j=1:3
                func(pt) = local_basis_func(C[:, i]', pt)*local_basis_func(C[:, j]', pt)
                M_elem[i, j] = gaussian_quad2_tri(func, p[t[k, :], :])
            end
        end

        K_elem = zeros(3, 3)
        for i=1:3
            for j=1:3
                func(pt) = C[3, i]*C[3, j]
                K_elem[i, j] = gaussian_quad2_tri(func, p[t[k, :], :])
            end
        end

        G_elem_bot = zeros(3, 3)
        edge_pts = [] # how many edge points does the element have
        for i=1:3
            if t[k, i] in ebot
                push!(edge_pts, i)
            end
        end
        if size(edge_pts, 1) == 2 # if two points lie on the bottom edge, integral is nonzero
            p1 = p[t[k, edge_pts[1]], :]' 
            p2 = p[t[k, edge_pts[2]], :]'
            pt(t) = (1 - t)*p1 + t*p2 # integrate along bottom edge line
            for i=1:3
                for j=1:3
                    if i != edge_pts[1] && i != edge_pts[2]
                        continue # point not on edge kills integral (true even without if statement)
                    else
                        func(t) = local_basis_func(C[:, i]', pt(t))*C[3, j]
                        G_elem_bot[i, j] = gaussian_quad2(func, 0, 1)
                    end
                end
            end
        end

        G_elem_top = zeros(3, 3)
        edge_pts = [] # how many edge points does the element have
        for i=1:3
            if t[k, i] in etop
                push!(edge_pts, i)
            end
        end
        if size(edge_pts, 1) == 2 # if two points lie on the bottom edge, integral is nonzero
            p1 = p[t[k, edge_pts[1]], :]' 
            p2 = p[t[k, edge_pts[2]], :]'
            point(t) = (1 - t)*p1 + t*p2 # integrate along bottom edge line
            for i=1:3
                for j=1:3
                    if i != edge_pts[1] && i != edge_pts[2]
                        continue # point not on edge kills integral (true even without if statement)
                    else
                        func(t) = local_basis_func(C[:, i]', point(t))*C[3, j]
                        G_elem_top[i, j] = gaussian_quad2(func, 0, 1)
                    end
                end
            end
        end

        rhs_elem = zeros(3)
        for i=1:3
            func(pt) = local_basis_func(C[:, i]', pt)*bx(pt)
            rhs_elem[i] = gaussian_quad2_tri(func, p[t[k, :], :])
        end

		## add to global system
		for i=1:3
			for j=1:3
                # first eqtn
                push!(A, (umap[1, t[k, i]], umap[1, t[k, j]],  M_elem[i, j]))
                push!(A, (umap[1, t[k, i]], umap[2, t[k, j]],  K_elem[i, j]))
                push!(A, (umap[1, t[k, i]], umap[2, t[k, j]], -G_elem_top[i, j]))
                # second eqtn
                push!(A, (umap[2, t[k, i]], umap[1, t[k, j]],  nu*G_elem_top[i, j]))
                push!(A, (umap[2, t[k, i]], umap[1, t[k, j]], -nu*G_elem_bot[i, j]))
                push!(A, (umap[2, t[k, i]], umap[1, t[k, j]], -nu*K_elem[i, j]))
                push!(A, (umap[2, t[k, i]], umap[2, t[k, j]],  nu*M_elem[i, j]))
			end
		end
        rhs[umap[2, t[k, :]]] -= rhs_elem
	end
    
    A = sparse((x->x[1]).(A), (x->x[2]).(A), (x->x[3]).(A), nPts*nVars, nPts*nVars)
  
    # bottom dirichlet
    for iVar=2
        for i=1:size(ebot, 1)
            swapEq = 2
            row = umap[swapEq, ebot[i]]
            A[row, :] .= 0
            A[row, umap[iVar, ebot[i]]] = 1
            rhs[row] = 0
        end
    end

    # top dirichlet
    for iVar=1:2
        for i=1:size(etop, 1)
            row = umap[iVar, etop[i]]
            A[row, :] .= 0
            A[row, row] = 1
            rhs[row] = 0
        end
    end

    # periodic bdy conditions
    for iVar=1:2
        for i=1:size(eleft, 1)
            rowLeft = umap[iVar, eleft[i]]
            rowRight = umap[iVar, eright[i]]

            A[rowLeft, :] += A[rowRight, :]
            rhs[rowLeft] += rhs[rowRight]

            A[rowRight, :] .= 0
            A[rowRight, rowRight] = 1
            A[rowRight, rowLeft] = -1
            rhs[rowRight] = 0
        end
    end

	dropzeros!(A)

    println("rank(A)    = ", rank(A))
    println("nPts*nVars = ", nPts*nVars)

	# Solve system
	sol = A \ rhs 
	sol = reshape(sol, (nVars, nPts))
	return sol
end

"""
solve problem and plot
"""
function solve(p, t, e)
    x = 0:L/1000:L
    topo = @. -H(x)

    nPts = size(p, 1)

    # buoyancy dist
    #= b(pt) = N^2*amp*exp(-(pt[2] + H(pt[1]))^2/2/(H0/4)^2) =#
    b(pt) = N^2*1e1*sin(2*pt[1]*pi/L)*exp(-(pt[2] + H0)^2/2/(H0/4)^2)
    #= b(pt) = N^2*1e1*exp(-(pt[2] + H0/2)^2/2/(H0/8)^2 - (pt[1] - L/2)^2/2/(L/8)^2) =#
    
    # x gradient of buoyancy dist
    #= bx(pt) = Hx(pt[1])*(pt[2] + H(pt[1]))/(H0/4)^4*b(pt) =#
    bx(pt) = N^2*1e1*2*pi/L*cos(2*pt[1]*pi/L)*exp(-(pt[2] + H0)^2/2/(H0/4)^2)
    #= bx(pt) = (pt[1] - L/2)/(L/8)^2*N^2*1e1*exp(-(pt[2] + H0/2)^2/2/(H0/8)^2 - (pt[1] - L/2)^2/2/(L/8)^2) =#

    rho = zeros(nPts)
    for i=1:nPts
        rho[i] = N^2*p[i, 2] + b(p[i, :])
    end
    
    # solve
	sol = fem_inversion(p, t, e, bx)
    omega = sol[1, :]
    chi = sol[2, :]

    # plot
    figure()
    vmax = maximum(abs.(omega))
    vmin = -vmax
    tplot(p, t, omega, cmap="RdBu_r", vmin=vmin, vmax=vmax)
    colorbar(label=L"$\omega$ (s$^{-1}$)")
    tricontour(p[:, 1], p[:, 2], t .- 1, rho, 10, colors="k", alpha=0.5)
    plot(x, topo, "k-")
    title(L"$\omega$")
    xlabel(L"$x$ (m)")
    ylabel(L"$z$ (m)")

    # plot
    figure()
    vmax = maximum(abs.(chi))
    vmin = -vmax
    tplot(p, t, chi, cmap="RdBu_r", vmin=vmin, vmax=vmax)
    colorbar(label=L"$\chi$ (m$^2$ s$^{-1}$)")
    tricontour(p[:, 1], p[:, 2], t .- 1, rho, 10, colors="k", alpha=0.5)
    plot(x, topo, "k-")
    title(L"$\chi$")
    xlabel(L"$x$ (m)")
    ylabel(L"$z$ (m)")
end

# generate mesh and load it as p, t, e data structure
lc = 2^-6
bdyRef = 2^0
#= generateMesh(lc, bdyRef) =#
#= error() =#
p, t, e = loadMesh()
#= tplot(p, t) =#

# find left, right, top, and bottom edge nodes
eleft = e[@. p[e, 1] == 0]
eleft = eleft[end:-1:1]
eleft = [eleft[end-1]; eleft[end]; eleft[1:end-2]]
eright = e[@. p[e, 1] == L]
etop = e[@. p[e, 2] == 0]
#= ebot = e[@. abs(p[e, 2] + H(p[e, 1]) < H0/100)] =#
ebot = e[@. (p[e, 1] != 0) & (p[e, 1] != L) & (p[e, 2] != 0)]
push!(ebot, e[@. (p[e, 1] == 0) & (p[e, 2] == -H(0))][1])
push!(ebot, e[@. (p[e, 1] == L) & (abs(p[e, 2] + H(L)) <= lc/2)][1])

# solve!
solve(p, t, e)
