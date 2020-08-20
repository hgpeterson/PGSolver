using SparseArrays, PyPlot, LinearAlgebra
include("setParams.jl")
include("mesh.jl")

#= close("all") =#
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
    return cs*[1 pt]'
end


"""
	u = fem_ekman(p, t, e)
Solves Ekman layer equations on domain described by unstructured
triangular mesh `p, t`. The boundary conditions are homogeneous Neumann except
for the nodes in the array `e` which are homogeneous Dirichlet.
"""
function fem_ekman(p, t, e)
	N = size(p, 1)
	K = size(t, 1)

    nVars = 2

    umap = reshape(1:(N*nVars), nVars, N) 

	# Create global linear system using stamping method
    A = Tuple{Int64,Int64,Float64}[]  
    b = zeros(N*nVars)
	for k = 1:K
		# Calculate Elementary matrices and loads

		## get coeffs for linear basis func c1 + c2 x + c3 z
		A_elem = zeros(3, 3)
		b_elem = zeros(3, 1)
		V = zeros(3, 3)
		for row = 1:3
			V[row, :] = [1 p[t[k, row], 1] p[t[k, row], 2]]
		end
		C = inv(V)

		## calculate entries of matrix
        M_elem = zeros(3, 3)
        for i=1:3
            for j=1:3
                func(pt) = f*local_basis_func(C[:, i]', pt)*local_basis_func(C[:, j]', pt)
                M_elem[i, j] = gaussian_quad2_tri(func, p[t[k, :], :])
            end
        end

        K_elem_z_nu = zeros(3, 3)
        for i=1:3
            for j=1:3
                func(pt) = nu(pt[1], pt[2])*C[3, i]*C[3, j]
                K_elem_z_nu[i, j] = gaussian_quad2_tri(func, p[t[k, :], :])
            end
        end

        G_elem_z_nu = zeros(3, 3)
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
                        G_elem_z_nu[i, j] = 0 # point not on edge kills integral (true even without if statement)
                    else
                        func(t) = nu(pt(t)[1], pt(t)[2])*local_basis_func(C[:, i]', pt(t))*C[3, j]
                        G_elem_z_nu[i, j] = gaussian_quad2(func, 0, 1)
                    end
                end
            end
        end

		## add to global system
		for i=1:3
			for j=1:3
                # first eqtn
                push!(A, (umap[1, t[k, i]], umap[2, t[k, j]], -M_elem[i, j]))
                push!(A, (umap[1, t[k, i]], umap[1, t[k, j]], K_elem_z_nu[i, j]))
                #= push!(A, (umap[1, t[k, i]], umap[1, t[k, j]], G_elem_z_nu[i, j])) =#
                # second eqtn
                push!(A, (umap[2, t[k, i]], umap[1, t[k, j]], M_elem[i, j]))
                push!(A, (umap[2, t[k, i]], umap[2, t[k, j]], K_elem_z_nu[i, j]))
                #= push!(A, (umap[2, t[k, i]], umap[2, t[k, j]], G_elem_z_nu[i, j])) =#
			end
		end
	end
    
    A = sparse((x->x[1]).(A), (x->x[2]).(A), (x->x[3]).(A), N*nVars, N*nVars)

    # bottom dirichlet
    for iVar=1:2
        for i=1:size(ebot, 1)
            row = umap[iVar, ebot[i]]
            A[row, :] .= 0
            A[row, row] = 1
            b[row] = 1e-2
        end
    end

    # right bdy is copy of left
    for iVar=1:2
        for i=1:size(eleft, 1)
            rowLeft = umap[iVar, eleft[i]]
            rowRight = umap[iVar, eright[i]]

            A[rowLeft, :] += A[rowRight, :]
            b[rowLeft] += b[rowRight]

            A[rowRight, :] .= 0
            A[rowRight, rowRight] = 1
            A[rowRight, rowLeft] = -1
            b[rowRight] = 0
        end
    end

	dropzeros!(A)

    println(rank(A))
    println(N*nVars)

	# Solve system
	sol = A \ b 
	sol = reshape(sol, (nVars, N))
	return sol
end

"""
solve problem and plot
"""
function solve(p, t, e)
    # solve
	sol = fem_ekman(p, t, e)
    u = sol[1, :]
    v = sol[2, :]

    # plot
    figure()
    umax = maximum(abs.(u))
    umin = -umax
    tplot(p, t, u, cmap="bwr", vmin=umin, vmax=umax)
    colorbar(label=L"$u$ (m s$^{-1}$)")
    title(L"$u$")
    xlabel(L"$y$ (m)")
    ylabel(L"$z$ (m)")

    figure()
    vmax = maximum(abs.(v))
    vmin = -vmax
    tplot(p, t, v, cmap="bwr", vmin=vmin, vmax=vmax)
    colorbar(label=L"$v$ (m s$^{-1}$)")
    title(L"$v$")
    xlabel(L"$y$ (m)")
    ylabel(L"$z$ (m)")
end

# generate mesh and load it as p, t, e data structure
lc = 2^-7
bdyRef = 2^3
generateMesh(lc, bdyRef)
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
