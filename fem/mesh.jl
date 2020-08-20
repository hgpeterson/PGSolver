include("meshutils.jl")
#= include("setParams.jl") =#
import gmsh

"""
    generateMesh(lc)

Creates 2D mesh of hill with characteristic length `lc`.
Data staved to "mesh.msh".
"""
function generateMesh(lc, bdyRef)
    # init
    gmsh.initialize()
    
    # log
    gmsh.option.setNumber("General.Terminal", 1)
    
    # model
    gmsh.model.add("hillMesh")
    
    # edge points
    #= bdyRef = 10 # boundary resolution factor =# 
    #= bdyRef = 1 =#
    pts = []
    push!(pts, gmsh.model.geo.addPoint(0, 0, 0, bdyRef*lc))
    push!(pts, gmsh.model.geo.addPoint(0, -H(0)/H0, 0, lc))
    for i=1:Int64(1/lc)
        x = i*lc
        push!(pts, gmsh.model.geo.addPoint(x, -H(x*L)/H0, 0, lc))
    end
    push!(pts, gmsh.model.geo.addPoint(1, 0, 0, bdyRef*lc))
    
    # connect edge points by lines
    curves = []
    for i=1:size(pts, 1)-1
        push!(curves, gmsh.model.geo.addLine(pts[i], pts[i+1]))
    end
    push!(curves, gmsh.model.geo.addLine(pts[end], pts[1]))
    
    # loop curves together and define surface
    gmsh.model.geo.addCurveLoop(curves, 1)
    gmsh.model.geo.addPlaneSurface([1], 1)
    
    # sync
    gmsh.model.geo.synchronize()
    
    # make left bdy copy of right
    translation = [1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
    gmsh.model.mesh.setPeriodic(1, [curves[1]], [curves[end-1]], translation)
    
    # sync and generate
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)
    
    # save and show
    gmsh.write("mesh.msh")
    gmsh.finalize()
end

"""
    p, t, e = loadMesh()

Loads mesh from "mesh.msh" file and returns data structure as in "pmesh.jl"
where `p` is an nNodes x nDims array of node positions, `t` is an
nTri x 3 array of triangle nodes indeices, and `e` is an nEdgeNodes x 1 array of
edge node indices.
"""
function loadMesh()
    # initialize mesh and load from file
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("hillMesh")
    gmsh.open("mesh.msh")

    # find triangle nodes from the elements in the type-2 surface with tag 1
    tri_nodes = gmsh.model.mesh.getElements(2, 1)[3][1]
    nTri = Int64(size(tri_nodes, 1)/3)
    t = zeros(nTri, 3)
    for i=1:nTri
        t[i, :] = [tri_nodes[3*i-2] tri_nodes[3*i-1] tri_nodes[3*i]]
    end
    t = Int64.(t)

    # find node positions by looping through indices
    nPts = Int64(maximum(t))
    p = zeros(nPts, 2)
    for i=1:nPts
        p[i, :] = gmsh.model.mesh.getNode(i)[1][1:2]
    end

    # rescale
    p[:, 1] .*= L
    p[:, 2] .*= H0

    # get edge nodes (see meshutils.jl)
    e = boundary_nodes(t)

    gmsh.finalize()

    return p, t, e
end
