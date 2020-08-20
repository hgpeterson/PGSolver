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

#= """ =#
#= compute elementary matrix =#
#= """ =#
#= function get_elem_mat(func, pts, area) =#
#=     A_elem = zeros(3, 3) =#
#=     for i=1:3 =#
#=         for j=1:3 =#
#=             A_elem[i, j] = area*gaussian_quad2(func, pts) =#
#=         end =#
#=     end =#
#= end =#

"""
	sol = fem_inversion(p, t, e, b)
Inverts buoyancy field `b` for flow field on domain described by unstructured
triangular mesh `p, t`. 
"""
function fem_inversion(p, t, e, b)
	nPts = size(p, 1)
	nTri = size(t, 1)

    nVars = 4

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

        #= # vectorized =#
        #= # integration pts and weights for gaussian quadrature =#
		#= integration_pts = ones(3, 3) =#
		#= integration_pts[2:3, :] = (1/6*sum(p[t[k, :], :], dims=1) .+ 1/2*I*p[t[k, :], :])' =#
        #= weight = 1/3 =#
		#= M_elem = f*area*weight*(C'*integration_pts)*(C'*integration_pts)' =#
        #= C_x = zeros(3, 3) =#
        #= C_x[1, :] = C[2, :] =#
		#= C_elem_x = area*weight*(C'*integration_pts)*(C_x'*integration_pts)' =#
        #= C_z = zeros(3, 3) =#
        #= C_z[1, :] = C[3, :] =#
		#= C_elem_z = area*weight*(C'*integration_pts)*(C_z'*integration_pts)' =#
		#= K_elem_z_kap = kap*area*weight*(C_z'*integration_pts)*(C_z'*integration_pts)' =#

        # for loops
        M_elem = zeros(3, 3)
        for i=1:3
            for j=1:3
                func(pt) = f*local_basis_func(C[:, i]', pt)*local_basis_func(C[:, j]', pt)
                M_elem[i, j] = gaussian_quad2_tri(func, p[t[k, :], :])
            end
        end

        C_elem_x = zeros(3, 3)
        for i=1:3
            for j=1:3
                func(pt) = local_basis_func(C[:, i]', pt)*C[2, j]
                C_elem_x[i, j] = gaussian_quad2_tri(func, p[t[k, :], :])
            end
        end

        C_elem_z = zeros(3, 3)
        for i=1:3
            for j=1:3
                func(pt) = local_basis_func(C[:, i]', pt)*C[3, j]
                C_elem_z[i, j] = gaussian_quad2_tri(func, p[t[k, :], :])
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

        rhs_elem = zeros(3)
        for i=1:3
            func(pt) = local_basis_func(C[:, i]', pt)*b(pt)
            rhs_elem[i] = gaussian_quad2_tri(func, p[t[k, :], :])
        end

		## add to global system
		for i=1:3
			for j=1:3
                # first eqtn
                push!(A, (umap[1, t[k, i]], umap[2, t[k, j]], -M_elem[i, j]))
                push!(A, (umap[1, t[k, i]], umap[4, t[k, j]], C_elem_x[i, j]))
                push!(A, (umap[1, t[k, i]], umap[1, t[k, j]], K_elem_z_nu[i, j]))
                push!(A, (umap[1, t[k, i]], umap[1, t[k, j]], G_elem_z_nu[i, j]))
                # second eqtn
                push!(A, (umap[2, t[k, i]], umap[1, t[k, j]], M_elem[i, j]))
                push!(A, (umap[2, t[k, i]], umap[2, t[k, j]], K_elem_z_nu[i, j]))
                push!(A, (umap[2, t[k, i]], umap[2, t[k, j]], G_elem_z_nu[i, j]))
                #= # third eqtn =#
                #= push!(A, (umap[3, t[k, i]], umap[4, t[k, j]], C_elem_z[i, j])) =#
                #= # fourth eqtn =#
                #= push!(A, (umap[4, t[k, i]], umap[1, t[k, j]], C_elem_x[i, j])) =#
                #= push!(A, (umap[4, t[k, i]], umap[3, t[k, j]], C_elem_z[i, j])) =#
                # third eqtn
                push!(A, (umap[3, t[k, i]], umap[1, t[k, j]], C_elem_x[i, j]))
                push!(A, (umap[3, t[k, i]], umap[3, t[k, j]], C_elem_z[i, j]))
                # fourth eqtn
                push!(A, (umap[4, t[k, i]], umap[4, t[k, j]], C_elem_z[i, j]))
			end
		end
        #= rhs[umap[3, t[k, :]]] .+= rhs_elem =#
        rhs[umap[4, t[k, :]]] .+= rhs_elem
	end
    
    A = sparse((x->x[1]).(A), (x->x[2]).(A), (x->x[3]).(A), nPts*nVars, nPts*nVars)
  
    # bottom dirichlet
    for iVar=1:3
    #= for iVar=1:4 =#
        for i=1:size(ebot, 1)
            row = umap[iVar, ebot[i]]
            A[row, :] .= 0
            A[row, row] = 1
            rhs[row] = 0
        end
    end

    # top dirichlet
    for i=1:size(etop, 1)
        #= row = umap[3, etop[i]] =#
        #= A[row, :] .= 0 =#
        #= A[row, row] = 1 =#
        #= rhs[row] = 0 =#
        eqtn = 3
        A[umap[eqtn, etop[i]], :] .= 0
        A[umap[eqtn, etop[i]], umap[3, etop[i]]] = 1
        rhs[umap[eqtn, etop[i]]] = 0
    end

    #= # make p zero at a point =#
    #= i = 1 =#
    #= row = umap[4, ebot[i]] =#
    #= A[row, :] .= 0 =#
    #= A[row, row] = 1 =#
    #= rhs[row] = 0 =#
    # make sum of p zero (replace eqtn 4 at a point)
    i = 10
    row = umap[4, ebot[i]]
    A[row, :] .= 0
    for j=1:nPts
        A[row, umap[4, j]] = 1
    end
    rhs[row] = 0
    #= tplot(p, t) =#
    #= plot(p[ebot[i], 1], p[ebot[i], 2], "ro") =#

    # periodic bdy conditions
    for iVar=1:4
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
    #= b(pt) = 0 =#
    #= b(pt) = N^2/100*(H(pt[1]) - H0) =#
    b(pt) = N^2*1e2*sin(2*pt[1]*pi/L)
    #= b(pt) = N^2*amp*exp(-(pt[2] + H(pt[1]))^2/2/(H0/4)^2) =#

    rho = zeros(nPts)
    for i=1:nPts
        rho[i] = N^2*p[i, 2] + b(p[i, :])
    end
    
    # solve
	sol = fem_inversion(p, t, e, b)
    u = sol[1, :]
    v = sol[2, :]
    w = sol[3, :]
    pres = sol[4, :]

    # plot
    figure()
    vmax = maximum(abs.(u))
    vmin = -vmax
    tplot(p, t, u, cmap="bwr", vmin=vmin, vmax=vmax)
    colorbar(label=L"$u$ (m s$^{-1}$)")
    tricontour(p[:, 1], p[:, 2], t .- 1, rho, 10, colors="k", alpha=0.5)
    plot(x, topo, "k-")
    title(L"$u$")
    xlabel(L"$x$ (m)")
    ylabel(L"$z$ (m)")

    figure()
    vmax = maximum(abs.(v))
    vmin = -vmax
    tplot(p, t, v, cmap="bwr", vmin=vmin, vmax=vmax)
    colorbar(label=L"$v$ (m s$^{-1}$)")
    tricontour(p[:, 1], p[:, 2], t .- 1, rho, 10, colors="k", alpha=0.5)
    plot(x, topo, "k-")
    title(L"$v$")
    xlabel(L"$x$ (m)")
    ylabel(L"$z$ (m)")

    figure()
    vmax = maximum(abs.(w))
    vmin = -vmax
    tplot(p, t, w, cmap="bwr", vmin=vmin, vmax=vmax)
    colorbar(label=L"$w$ (m s$^{-1}$)")
    tricontour(p[:, 1], p[:, 2], t .- 1, rho, 10, colors="k", alpha=0.5)
    plot(x, topo, "k-")
    title(L"$w$")
    xlabel(L"$x$ (m)")
    ylabel(L"$z$ (m)")

    figure()
    tplot(p, t, pres)
    colorbar(label=L"$p$ (m$^2$ s$^{-2}$)")
    tricontour(p[:, 1], p[:, 2], t .- 1, rho, 10, colors="k", alpha=0.5)
    plot(x, topo, "k-")
    title(L"$p$")
    xlabel(L"$x$ (m)")
    ylabel(L"$z$ (m)")
end

# generate mesh and load it as p, t, e data structure
lc = 2^-7
bdyRef = 2^3
generateMesh(lc, bdyRef)
#= error() =#
p, t, e = loadMesh()
tplot(p, t)
error()

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
