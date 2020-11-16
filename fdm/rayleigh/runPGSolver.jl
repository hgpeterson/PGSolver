using PyPlot, PyCall, Printf, SparseArrays, LinearAlgebra, HDF5, SpecialFunctions

plt.style.use("~/paper_plots.mplstyle")
close("all")
pygui(false)

include("../../myJuliaLib.jl")
include("setParams.jl")
include("terrainFollowing.jl")
include("plottingLib.jl")
include("inversion.jl")
include("evolution.jl")

################################################################################
# run evolution integrations
################################################################################

#= # integrate =#
#= b = evolve(500) =#
#= profilePlot(["b1000.h5", "b2000.h5", "b3000.h5", "b4000.h5", "b5000.h5"], 1) =#


################################################################################
# make some plots
################################################################################

#= # 1D sol =#
#= for t=[1000 2000 3000 4000 5000 Inf] =#
#=     b1D, u1D, v1D, w1D = pointwise1D(t*86400) =#
#=     ridgePlot(b1D, b1D, "", "") =#
#=     savefig(string("b1D", t, ".png")) =#
#= end =#

#= # redo profile plots =#
#= folder = "images/constKappa/full2D/" =#
#= profilePlot(string.(folder, ["b1000.h5", "b2000.h5", "b3000.h5", "b4000.h5", "b5000.h5"]), 1) =#

#= # bξ vs bσ =#
#= for bfile in string.("images/constKappa/full2D/", ["b1000.h5", "b2000.h5", "b3000.h5", "b4000.h5", "b5000.h5"]) =#
#=     file = h5open(bfile, "r") =#
#=     b = read(file, "b") =#
#=     t = read(file, "t") =#
#=     close(file) =#

#=     ridgePlot(ξDerivativeTF(b), b, @sprintf("t = %d days", t/86400), L"$b_\xi$") =#
#=     savefig(@sprintf("b_xi%d.png", t/86400)) =#
#=     ridgePlot(σσ.*Hx.(ξξ)./H.(ξξ).*σDerivativeTF(b), b, @sprintf("t = %d days", t/86400), L"$\sigma H_xb_\sigma/H$") =#
#=     savefig(@sprintf("b_sig%d.png", t/86400)) =#
#= end =#

# px
inversionLHS = lu(getInversionLHS())
folder = "images/constKappa/full2D/"
for bfile in string.(folder, ["b1000.h5", "b2000.h5", "b3000.h5", "b4000.h5", "b5000.h5"])
    file = h5open(bfile, "r")
    b = read(file, "b")
    t = read(file, "t")
    close(file)

    chi, uξ, uη, uσ, U = invert(b, inversionLHS) 

    u, v, w = transformFromTF(uξ, uη, uσ)

    px = f*v + r*u

    ridgePlot(px, b, @sprintf("t = %d days", t/86400), L"$p_x$")
    savefig(@sprintf("p_x%d.png", t/86400))
end
