include("setParams.jl")
include("myJuliaLib.jl")
include("plottingLib.jl")
include("inversion.jl")
include("evolution.jl")

################################################################################
# manually set `b` if you want
################################################################################
#= b = zeros(nξ, nσ) =#
#= bx = zeros(nξ, nσ) =#
#= for i=1:nξ =#
#=     for j=1:nσ =#
#=         # exponential in σ (centered at bottom) =#
#=         decay_scale = 200/H0 =#
#=         A = decay_scale*N^2*H0 =#
#=         b[i, j] = A*exp(-N^2*H(ξ[i])*(σ[j] + 1)/A) =#
#=         bx[i, j] = -N^2*Hx(ξ[i])*exp(-N^2*H(ξ[i])*(σ[j] + 1)/A) =#

#=         #1= # gaussian in σ (centered at bottom) =1# =#
#=         #1= b[i, j] = N^2*amp*exp(-(σ[j] + 1)^2/2/(1/4)^2) =1# =#
#=         #1= bx[i, j] = N^2*amp*exp(-(σ[j] + 1)^2/2/(1/4)^2)*(σ[j] + 1)/(1/4)^2*σ[j]*Hx(ξ[i])/H(ξ[i]) =1# =#
#=     end =#
#= end =#

#= b = @. h*N^2*exp(-(z + H(x))/h) =#
#= #1= bx = @. -N^2*Hx(x)*exp(-(z + H(x))/h) =1# =#
#= bx = zeros(nξ, nσ) =#
#= for j=1:nσ =#
#=     bx[:, j] = differentiate(b[:, j], ξ) =#
#= end =#
#= for i=1:nξ =#
#=     bx[i, :] .-= Hx(ξ[i])*σ.*differentiate(b[i, :], σ)/H(ξ[i]) =#
#= end =#

################################################################################
# test 1DAdjusted inversion
################################################################################
#= inversionLHS = getInversionLHS() =#
#= chi, uξ, uη, uσ, U = invert(b, inversionLHS) =#
#= u, v, w = transformFromTF(uξ, uη, uσ) =#
#= ridgePlot(chi, b, "chi", "") =#
#= savefig("chi.png") =#
#= close() =#
#= ridgePlot(u, b, "u", "") =#
#= savefig("u.png") =#
#= close() =#
#= ridgePlot(v, b, "v", "") =#
#= savefig("v.png") =#
#= close() =#
#= ridgePlot(w, b, "w", "") =#
#= savefig("w.png") =#
#= close() =#

#= inversionLHS1DAdjusted = getInversionLHS1DAdjusted() =#
#= inversionRHS1DAdjusted = getInversionRHS1DAdjusted(b) =#
#= sol = inversionLHS1DAdjusted\inversionRHS1DAdjusted =#
#= chi1D, u1D, v1D, w1D = postProcess1DAdjusted(sol) =#
#= ridgePlot(chi1D, b, "chi 1D", ""; vext=maximum(abs.(chi))) =#
#= savefig("chi1D_U.png") =#
#= close() =#
#= ridgePlot(u1D, b, "u 1D", ""; vext=maximum(abs.(u))) =#
#= savefig("u1D_U.png") =#
#= close() =#
#= ridgePlot(v1D, b, "v 1D", ""; vext=maximum(abs.(v))) =#
#= savefig("v1D_U.png") =#
#= close() =#
#= ridgePlot(w1D, b, "w 1D", ""; vext=maximum(abs.(w))) =#
#= savefig("w1D_U.png") =#
#= close() =#

################################################################################
# run evolution integrations
################################################################################

#= b = evolve(500) =#
#= profilePlot(["b1000.h5", "b2000.h5", "b3000.h5", "b4000.h5", "b5000.h5"], 1) =#

#= b = evolve1DAdjusted(500) =#
profilePlot1DAdjusted(["b1000.h5", "b2000.h5", "b3000.h5", "b4000.h5", "b5000.h5"], 1)
