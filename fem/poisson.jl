using SparseArrays, PyPlot, LinearAlgebra
include("setParams.jl")
include("mesh.jl")

close("all")
plt.style.use("~/presentation_plots.mplstyle")

"""
comupte area of triangle defined by pts
"""
function tri_area(pts)
	p1 = pts[1, :]; p2 = pts[2, :]; p3 = pts[3, :]
	area = 1/2 * abs(p1[1]*(p2[2]-p3[2]) + p2[1]*(p3[2]-p1[2]) + p3[1]*(p1[2]-p2[2]))
    return area
end

#= """ =#
#= 	sol = fem_inversion(p, t, e, b) =#
#= Inverts buoyancy field `b` for flow field on domain described by unstructured =#
#= triangular mesh `p, t`. =# 
#= """ =#
function fem_inversion(p, t, e)
	nPts = size(p, 1)
	nTri = size(t, 1)

    nVars = 1
    f = 1

    umap = reshape(1:(nPts*nVars), nVars, nPts) 

	# Create global linear system using stamping method
    A = Tuple{Int64,Int64,Float64}[]  
    b = zeros(nPts*nVars)
	for k = 1:nTri
		# Calculate Elementary matrices and loads

		## get coeffs for linear basis func c1 + c2 x + c3 z
		V = zeros(3, 3)
		for row = 1:3
			V[row, :] = [1 p[t[k, row], 1] p[t[k, row], 2]]
		end
		C = inv(V)

		## calculate entries of matrix
        area = tri_area(p[t[k, :], :])
		A_elem = area * (C[2,:]*C[2,:]' + C[3,:]*C[3,:]')
		b_elem = area / 3 * f 

		## add to global system
		for i=1:3
			for j=1:3
                push!(A, (umap[1, t[k, i]], umap[1, t[k, j]], A_elem[i, j]))
			end
		end
        b[umap[1, t[k, :]]] .+= b_elem
	end
    
    A = sparse((x->x[1]).(A), (x->x[2]).(A), (x->x[3]).(A), nPts*nVars, nPts*nVars)

    tplot(p, t)
    # bottom dirichlet
    for iVar=1:1
        for i=1:size(ebot, 1)
            #= plot(p[ebot[i], 1], p[ebot[i], 2], "ro") =#
            row = umap[iVar, ebot[i]]
            A[row, :] .= 0
            A[row, row] = 1
            b[row] = 0
        end
    end

    # periodic bdy conditions on left and right
    for iVar=1:1
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
    println(nPts*nVars)

	# Solve system 
    sol = A \ b 
    sol = reshape(sol, (nVars, nPts))
	return sol
end

"""
solve problem and plot
"""
function solve(p, t, e)
    nPts = size(p, 1)

    # solve
	sol = fem_inversion(p, t, e)
    u = sol[1, :]
    figure()
    tplot(p, t, u)
end

# generate mesh and load it as p, t, e data structure
lc = 2^-7
bdyRef = 10
#= generateMesh(lc, bdyRef) =#
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
