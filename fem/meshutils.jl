# Various mesh utilities, mostly for unstructured triangular meshes
#
# UC Berkeley Math 228B, Per-Olof Persson <persson@berkeley.edu>

using PyPlot, PyCall

"""
    t = delaunay(p)

Delaunay triangulation `t` of N x 2 node array `p`.
"""
function delaunay(p)
    tri = pyimport("matplotlib.tri")
    t = tri.Triangulation(p[:,1], p[:,2])
    return Int64.(t.triangles .+ 1)
end

#= """ =#
#=     edges, boundary_indices = all_edges(t) =#

#= Find all unique edges in the triangulation `t` (ne x 2 array) =#
#= Second output is indices to the boundary edges. =#
#= """ =#
#= function all_edges(t) =#
#=     edges = vcat(t[:,[1,2]], t[:,[2,3]], t[:,[3,1]]) =#
#=     edges = sortslices(sort(edges, dims=2), dims=1) =#
#=     dup = all(edges[2:end,:] - edges[1:end-1,:] .== 0, dims=2)[:] =#
#=     keep = .![false;dup] =#
#=     edges = edges[keep,:] =#
#=     dup = [dup;false] =#
#=     dup = dup[keep] =#
#=     return edges, findall(.!dup) =#
#= end =#

"""
    edges, boundary_indices, emap = all_edges(t)

Find all unique edges in the triangulation `t` (ne x 2 array)
Second output is indices to the boundary edges.
Third output emap (nt x 3 array) is a mapping from local triangle edges
to the global edge list, i.e., emap[it,k] is the global edge number
for local edge k (1,2,3) in triangle it.
"""
function all_edges(t)
    etag = vcat(t[:,[1,2]], t[:,[2,3]], t[:,[3,1]])
    etag = hcat(sort(etag, dims=2), 1:3*size(t,1))
    etag = sortslices(etag, dims=1)
    dup = all(etag[2:end,1:2] - etag[1:end-1,1:2] .== 0, dims=2)[:]
    keep = .![false;dup]
    edges = etag[keep,1:2]
    emap = cumsum(keep)
    invpermute!(emap, etag[:,3])
    emap = reshape(emap,:,3)
    dup = [dup;false]
    dup = dup[keep]
    bndix = findall(.!dup)
    return edges, bndix, emap
end

"""
    e = boundary_nodes(t)

Find all boundary nodes in the triangulation `t`.
"""
function boundary_nodes(t)
    edges, boundary_indices = all_edges(t)
    return unique(edges[boundary_indices,:][:])
end

"""
    tplot(p, t, u=nothing)

If `u` == nothing: Plot triangular mesh with nodes `p` and triangles `t`.
If `u` == solution vector: Plot filled contour color plot of solution `u`.
"""
function tplot(p, t, u=nothing; cmap="viridis", vmin=nothing, vmax=nothing)
    #= clf() =#
    #= axis("equal") =#
    if u == nothing
        tripcolor(p[:,1], p[:,2], t .- 1, 0*t[:,1],
                  cmap="Set3", edgecolors="k", linewidth=1)
    else
        if vmin == nothing
            vmin = minimum(u)
        end
        if vmax == nothing
            vmax = maximum(u)
        end
        #= tricontourf(p[:,1], p[:,2], t .- 1, u, 20) =#
        #= tripcolor(p[:,1], p[:,2], t .- 1, u, cmap=cmap, vmin=vmin, vmax=vmax, shading="flat") =#
        tripcolor(p[:,1], p[:,2], t .- 1, u, cmap=cmap, vmin=vmin, vmax=vmax, shading="gouraud")
    end
    #= draw() =#
end

"""
    inside = inpolygon(p, pv)

Determine if each point in the N x 2 node array `p` is inside the polygon
described by the NE x 2 node array `pv`.
"""
function inpolygon(p::Array{Float64,2}, pv::Array{Float64,2})
    path = pyimport("matplotlib.path")
    poly = path.Path(pv)
    inside = [poly.contains_point(p[ip,:]) for ip = 1:size(p,1)]
end

#= function test_inpolygon() =#
#=     pv = Float64[0 0; 1 0; 1 1; 0 1; 0.5 0.6; 0 0] =#
#=     x = collect(-.15:0.1:1.15) =#
#=     p = [x 0*x.+0.5] =#
#=     inside = inpolygon(p, pv) =#
#=     clf() =#
#=     plot(pv[:,1], pv[:,2], =#
#=          p[inside,1], p[inside,2], "*", =#
#=          p[.!inside,1], p[.!inside,2], "o") =#
#=     return =#
#= end =#
