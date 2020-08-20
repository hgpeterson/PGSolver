using SparseArrays, PyPlot, LinearAlgebra
include("setParams.jl")
include("mesh.jl")

close("all")
plt.style.use("~/presentation_plots.mplstyle")

"""
	p2, t2, e2 = p2mesh(p, t)

New pmesh for quadratic functions which produces
`p2`: N x 2, node coords of original mesh plus new midpoints
`t2`: T x 6, a local-to-global mapping for the T triangle elements
`e2`: E x 1, indices of boundary nodes for original mesh plus new midpoints
"""
function p2mesh(p, t)
	# Find all the edges at first
    edges, boundary_indices, emap = all_edges(t)

	# Add the midpoints of each edge
    midpts = 1/2 * reshape(p[edges[:, 1], :] + p[edges[:, 2], :], (size(edges, 1), 2))
    p2 = [p; midpts]

	# Add the midpoints of each triangle
    t2 = zeros(size(t, 1), 6)
    t2[:, 1:3] = t
    t2[:, 4:6] = size(p, 1) .+ emap
    t2 = convert.(Int64, t2)
    
    # Add the midpoints that were on the boundary
    e2 = [unique(edges[boundary_indices, :][:]); size(p, 1) .+ boundary_indices]
	return p2, t2, e2
end

"""
compute area of triangle defined by pts
"""
function tri_area(pts)
	x1 = pts[1, :]; x2 = pts[2, :]; x3 = pts[3, :]
	area = 1/2 * abs(x1[1]*(x2[2] - x3[2]) + x2[1]*(x3[2] - x1[2]) + x3[1]*(x1[2] - x2[2]))
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
evaluate local basis function defined by cs = [c1 c2 c3 c4 c5 c6] at pt = [x z]
"""
function local_basis_func(cs, pt)
    value = cs*[1 pt[1] pt[2] pt[1]^2 pt[1]*pt[2] pt[2]^2]'
    return value[1]
end

"""
	sol = fem_inversion(p, t, e, b)
Inverts buoyancy field `b` for flow field on domain described by unstructured
triangular mesh `p, t`. 
"""
function fem_inversion(p2, t2, e2, bx)
	nPts = size(p2, 1)
	nTri = size(t2, 1)

    iStress = nPts + 1

	# Create global linear system using stamping method
    A = Tuple{Int64,Int64,Float64}[]  
    rhs = zeros(nPts + 1)
	for k = 1:nTri
		# Calculate Elementary matrices and loads

		## get coeffs for quadratic basis func 
        # points in the element on the bottom boundary
        bot_pts = [] 
        for i=1:6
            if t2[k, i] in ebot
                push!(bot_pts, i)
            end
        end
        # van der Monde matrix (V*C = I)
		V = ones(6, 6)
		for row = 1:6
            x = p2[t2[k, row], 1]
            z = p2[t2[k, row], 2]
            # phi(x, z) = c1 + c2 x + c3 z + c4 x^2 + c5 xz + c6 z^2
			V[row, 2:end] = [x z x^2 x*z z^2]
		end
        #= RHS = Diagonal(ones(6)) =#
        # dz(phi) = 0 at bottom instead of phi = 1
        #= for i in bot_pts =#
        #=     #1= x = V[i, 2] =1# =#
        #=     #1= z = V[i, 3] =1# =#
        #=     #1= # dz(phi) = c3 + c5 x + 2 c6 z =1# =#
        #=     #1= V[i, :] = [0 0 1 0 x 2*z] =1# =#
        #=     # dz(phi) = 0 =#
        #=     RHS[i, :] .= 0 =#
        #= end =#
        #= display(V) =#
        #= display(RHS) =#
		#= C = V\RHS =#
        C = inv(V)

		## calculate entries of matrix
        K_elem = zeros(6, 6)
        for i=1:6
            for j=1:6
                func(pt) = nu*(2*C[6, i])*(2*C[6, j])
                K_elem[i, j] = gaussian_quad2_tri(func, p2[t2[k, 1:3], :])
            end
        end

        G_elem = zeros(6, 6)
        edge_pts = [] # how many edge points does the element have
        for i=1:6
            if t2[k, i] in ebot
                push!(edge_pts, i)
            end
        end
        if size(edge_pts, 1) == 3 # if three points lie on the bottom edge, integral is nonzero
            for i=1:6
                for j=1:6
                    if i in edge_pts && j in edge_pts
                        pt1 = p2[t2[k, i], :]' 
                        pt2 = p2[t2[k, j], :]'
                        pt(t) = (1 - t)*pt1 + t*pt2 # integrate along bottom edge line
                        func(t) = nu*(C[3, i] + C[5, i]*pt(t)[1] + 2*C[6, i]*pt(t)[2])*(2*C[6, j])
                        G_elem[i, j] = gaussian_quad2(func, 0, 1)
                    end
                end
            end
            #= display(G_elem) =#
        end

        M_elem = zeros(6, 6)
        for i=1:6
            for j=1:6
                func(pt) = f^2/nu*local_basis_func(C[:, i]', pt)*local_basis_func(C[:, j]', pt)
                M_elem[i, j] = gaussian_quad2_tri(func, p2[t2[k, 1:3], :])
            end
        end

        rhs_elem = zeros(6)
        for i=1:6
            func(pt) = bx(pt)*local_basis_func(C[:, i]', pt)
            rhs_elem[i] = gaussian_quad2_tri(func, p2[t2[k, 1:3], :])
        end

		## add to global system
		for i=1:6
			for j=1:6
                push!(A, (t2[k, i], t2[k, j], K_elem[i, j] + G_elem[i, j] + M_elem[i, j]))
			end
		end
        rhs[t2[k, :]] -= rhs_elem
	end
    
    # bottom stress
    push!(A, (iStress, iStress, 1.0))
    #= rhs[iStress] = 10 =#

    # space matrix
    A = sparse((x->x[1]).(A), (x->x[2]).(A), (x->x[3]).(A), nPts + 1, nPts + 1)
  
    # bottom dirichlet
    for i=1:size(ebot, 1)
        row = ebot[i]
        A[row, :] .= 0
        A[row, row] = 1
        rhs[row] = 0
    end

    # top dirichlet
    for i=1:size(etop, 1)
        row = etop[i]
        A[row, :] .= 0
        A[row, row] = 1
        A[row, iStress] = -1 # chi(0) = s
        rhs[row] = 0
    end

    # periodic 
    for i=1:size(eleft, 1)
        rowLeft = eleft[i]
        rowRight = eright[i]

        A[rowLeft, :] += A[rowRight, :]
        rhs[rowLeft] += rhs[rowRight]

        A[rowRight, :] .= 0
        A[rowRight, rowRight] = 1
        A[rowRight, rowLeft] = -1
        rhs[rowRight] = 0
    end
    
	dropzeros!(A)
    println("Matrix assembled.")

    #= println("rank(A)  = ", rank(A)) =#
    #= println("nPts + 1 = ", nPts + 1) =#

	# Solve system
	sol = A \ rhs 
	return sol
end

"""
solve problem and plot
"""
function solve(p2, t2, e2)
    x = 0:L/1000:L
    topo = @. -H(x)

    # buoyancy dist
    b(pt) = N^2*amp*exp(-(pt[2] + H(pt[1]))^2/2/(H0/4)^2)
    #= b(pt) = N^2*1e1*sin(2*pt[1]*pi/L)*exp(-(pt[2] + H0)^2/2/(H0/4)^2) =#
    #= b(pt) = N^2*1e1*exp(-(pt[2] + H0/2)^2/2/(H0/8)^2 - (pt[1] - L/2)^2/2/(L/8)^2) =#
    
    # x gradient of buoyancy dist
    bx(pt) = Hx(pt[1])*(pt[2] + H(pt[1]))/(H0/4)^4*b(pt)
    #= bx(pt) = N^2*1e1*2*pi/L*cos(2*pt[1]*pi/L)*exp(-(pt[2] + H0)^2/2/(H0/4)^2) =#
    #= bx(pt) = (pt[1] - L/2)/(L/8)^2*N^2*1e1*exp(-(pt[2] + H0/2)^2/2/(H0/8)^2 - (pt[1] - L/2)^2/2/(L/8)^2) =#

    rho = zeros(size(p, 1))
    for i=1:size(p, 1)
        rho[i] = N^2*p[i, 2] + b(p[i, :])
    end
    
    # solve
	sol = fem_inversion(p2, t2, e2, bx)
    chi2 = sol[1:end-1]
    chi = chi2[1:size(p, 1)]

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

    #= figure() =#
    #= vmax = maximum(abs.(v)) =#
    #= vmin = -vmax =#
    #= tplot(p, t, v, cmap="RdBu_r", vmin=vmin, vmax=vmax) =#
    #= colorbar(label=L"$v$ (m s$^{-1}$)") =#
    #= tricontour(p[:, 1], p[:, 2], t .- 1, rho, 10, colors="k", alpha=0.5) =#
    #= plot(x, topo, "k-") =#
    #= title(L"$v$") =#
    #= xlabel(L"$x$ (m)") =#
    #= ylabel(L"$z$ (m)") =#

    #= figure() =#
    #= vmax = maximum(abs.(w)) =#
    #= vmin = -vmax =#
    #= tplot(p, t, w, cmap="RdBu_r", vmin=vmin, vmax=vmax) =#
    #= colorbar(label=L"$w$ (m s$^{-1}$)") =#
    #= tricontour(p[:, 1], p[:, 2], t .- 1, rho, 10, colors="k", alpha=0.5) =#
    #= plot(x, topo, "k-") =#
    #= title(L"$w$") =#
    #= xlabel(L"$x$ (m)") =#
    #= ylabel(L"$z$ (m)") =#

    #= figure() =#
    #= tplot(p, t, pres) =#
    #= colorbar(label=L"$p$ (m$^2$ s$^{-2}$)") =#
    #= tricontour(p[:, 1], p[:, 2], t .- 1, rho, 10, colors="k", alpha=0.5) =#
    #= plot(x, topo, "k-") =#
    #= title(L"$p$") =#
    #= xlabel(L"$x$ (m)") =#
    #= ylabel(L"$z$ (m)") =#
end

# generate mesh and load it as p, t, e data structure
lc = 2^-5
bdyRef = 2^0
generateMesh(lc, bdyRef)
p, t, e = loadMesh()
#= tplot(p, t) =#

# add midpoints
p2, t2, e2 = p2mesh(p, t)

# find left, right, top, and bottom edge nodes
eleft = e2[@. p2[e2, 1] == 0]
eleft = eleft[end:-1:1]
eleft = [eleft[end-1]; eleft[end]; eleft[1:end-2]]
eright = e2[@. p2[e2, 1] == L]
etop = e2[@. p2[e2, 2] == 0]
ebot = e2[@. (p2[e2, 1] != 0) & (p2[e2, 1] != L) & (p2[e2, 2] != 0)]
push!(ebot, e2[@. (p2[e2, 1] == 0) & (p2[e2, 2] == -H(0))][1])
push!(ebot, e2[@. (p2[e2, 1] == L) & (abs(p2[e2, 2] + H(L)) <= lc/2)][1])

# solve!
solve(p2, t2, e2)
