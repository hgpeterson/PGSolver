# parameters

# canonical or transport-constrained case?
canonical = true
#= canonical = false =#

# RC20
H0 = 2e3
Pr = 1e0
f = -5.5e-5
N = 1e-3
θ = 2.5e-3
κ0 = 1e-4
κ1 = 1e-4
h = 100
v0 = 0.01

#= # RC20 - larger slope and v0 =#
#= H0 = 2e3 =#
#= Pr = 1e0 =#
#= f = -5.5e-5 =#
#= N = 1e-3 =#
#= θ = 1e-2 =#
#= κ0 = 1e-4 =#
#= κ1 = 1e-4 =#
#= h = 100 =#
#= v0 = 0.1 =#

#= # MR93 =#
#= H0 = 2e3 =#
#= Pr = 1e0 =#
#= f = 1e-4 =#
#= N = 3.5e-3 =#
#= θ = asin(1e-2) =#
#= κ0 = 1e-4 =#
#= κ1 = 1e-4 =#
#= h = 10 =#
#= v0 = -0.1 =#        

#= # MR91 =#
#= H0 = 1e0 =#
#= Pr = 1e3 =#
#= f = 1 =#
#= N = 2 =#
#= θ = 10*pi/180 =#
#= κ0 = 1e-9 =#
#= κ1 = 1e-9 =#
#= h = 10 =#
#= v0 = -0.01 =#        
#= Δt = 0.5 =#
#= tSave = 1 =#

# slope burger
S = N^2*tan(θ)^2/f^2

# canonical spin down timescale
τ = 1/abs(S*f*cos(θ))
Δt = τ/100
tSave = τ

# number of grid points
nẑ = 2^9

# grid
ẑ = @. H0*(1 - cos(pi*(0:nẑ-1)/(nẑ-1)))/2 # chebyshev (ẑ = 0 is bottom)

#= bottomIntense = true =#
bottomIntense = false
if bottomIntense
    κ = @. κ0 + κ1*exp(-(ẑ + H0)/h)
else
    κ = κ1*ones(nẑ)
end

# timestepping
adaptiveTimestep = false
α = 0.0

# Ekman layer
δ = sqrt(2*Pr*κ1/abs(f))

"""
    log(ofile, text)

Write `text` to `ofile` and print it.
"""
function log(ofile::IOStream, text::String)
    write(ofile, string(text, "\n"))
    println(text)
end

# log properties
ofile = open("out.txt", "w")
log(ofile, "\nSpin Down with Parameters\n")

log(ofile, @sprintf("nẑ = %1.5e", nẑ))
log(ofile, @sprintf("H0 = %1.5e m", H0))
log(ofile, @sprintf("Pr = %1.5e", Pr))
log(ofile, @sprintf("f  = %1.5e s-1", f))
log(ofile, @sprintf("N  = %1.5e s-1", N))
log(ofile, @sprintf("θ  = %1.5e rad", θ))
log(ofile, @sprintf("κ0 = %1.5e m2 s-1", κ0))
log(ofile, @sprintf("κ1 = %1.5e m2 s-1", κ1))
log(ofile, @sprintf("h  = %1.5e m", h))
log(ofile, @sprintf("v0 = %1.5e m s-1", v0))
log(ofile, @sprintf("Δt = %1.5e s = %1.5e days", Δt, Δt/86400))
log(ofile, @sprintf("α  = %1.5e", α))
log(ofile, @sprintf("S  = %1.5e", S))

log(ofile, string("Canonical:              ", canonical))
log(ofile, string("Bottom intensification: ", bottomIntense))
log(ofile, string("Adaptive timestep:      ", adaptiveTimestep))

log(ofile, @sprintf("\nEkman layer thickness ~ %1.5e m", δ))
log(ofile, @sprintf("          ẑ[2] - ẑ[1] ~ %1.5e m", ẑ[2] - ẑ[1]))

log(ofile, @sprintf("\nSpin-down ~ %1.5e s ~ %1.5e days\n", τ, τ/86400))
close(ofile)
