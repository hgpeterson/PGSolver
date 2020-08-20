using PyPlot
include("../meshutils.jl")
import gmsh

# params
L = 1e6
H0 = 1e3
H(x) = H0 - 4e2*sin(2*pi*x/L)

function save_mesh()
    # init
    gmsh.initialize()
    
    # log
    gmsh.option.setNumber("General.Terminal", 1)
    
    # model
    gmsh.model.add("t1")
    
    # default characteristic mesh length
    lc = 1e-2
    
    # edge points
    pts = []
    push!(pts, gmsh.model.geo.addPoint(0, 0, 0, lc))
    push!(pts, gmsh.model.geo.addPoint(0, -1, 0, lc))
    for i=1:Int64(1/lc)
        x = i*lc
        push!(pts, gmsh.model.geo.addPoint(x, -H(x*L)/H0, 0, lc))
    end
    push!(pts, gmsh.model.geo.addPoint(1, 0, 0, lc))
    
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
    gmsh.write("hill.msh")
    #= gmsh.fltk.run() =#
    gmsh.finalize()
end

function load_mesh()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("t1")
    gmsh.open("hill.msh")

    #= segments = gmsh.model.mesh.getElements(1,1)[2][1] =#
    #= println("Segment tags composing Line(1)") =#
    #= println(Int64.(segments)) =#
    #= nodes = gmsh.model.mesh.getElements(1,1)[3][1] =#
    #= nodes_per_elem = gmsh.model.mesh.getElementProperties(1)[4] =#
    #= println("Number of nodes per segment: ", nodes_per_elem) =#
    #= for j=1:size(segments, 1) =#
    #=     println("Element Tag(", segments[j], ")=[", nodes[2*j-1], ", ", nodes[2*j], "]") =#
    #= end =#

    #= segments = gmsh.model.mesh.getElements(2, 1)[2][1] =#
    #= println("Triangle tags composing Surface(1)") =#
    #= println(Int64.(segments)) =#
    #= nodes = gmsh.model.mesh.getElements(2, 1)[3][1] =#
    #= nodes_per_elem = gmsh.model.mesh.getElementProperties(2)[4] =#
    #= println("Number of nodes per triangle: ", nodes_per_elem) =#
    #= for j=1:size(segments, 1) =#
    #=     println("Element Tag(", segments[j], ")=[", nodes[3*j-2], ", ", nodes[3*j-1], ", ", nodes[3*j], "]") =#
    #= end =#

    tri_nodes = gmsh.model.mesh.getElements(2, 1)[3][1]
    nTri = Int64(size(tri_nodes, 1)/3)
    t = zeros(nTri, 3)
    for i=1:nTri
        t[i, :] = [tri_nodes[3*i-2] tri_nodes[3*i-1] tri_nodes[3*i]]
    end
    t = Int64.(t)

    nPts = Int64(maximum(t))
    p = zeros(nPts, 2)
    for i=1:nPts
        p[i, :] = gmsh.model.mesh.getNode(i)[1][1:2]
    end
    p[:, 1] .*= L
    p[:, 2] .*= H0

    e = boundary_indices(t)

    gmsh.finalize()

    return p, t, e
end

#= save_mesh() =#

p, t, e = load_mesh()
#= tplot(p, t) =#
