# parameters (as in RC20)
H0 = 2e3
Pr = 1e0

#= f = -5.5e-5 =#
#= N = 1e-3 =#
#= θ = 2.5e-3 =#
f = 1e-4
N = 3.5e-3
θ = asin(1e-2)

S = N^2*tan(θ)^2/f^2

# canonical or transport-constrained case?
canonical = true
#= canonical = false =#

# number of grid points
nẑ = 2^9

# grid
#= z = @. -H0*(cos(pi*(0:nẑ-1)/(nẑ-1)) + 1)/2 # chebyshev =# 
#= ẑ = z/cos(θ) =#
ẑ = @. -H0*(cos(pi*(0:nẑ-1)/(nẑ-1)) + 1)/2 # chebyshev 

# diffusivity
#= κ0 = 6e-5 =#
#= κ1 = 2e-3 =#
#= h = 200 =#
κ0 = 1e-4
κ1 = 1e-2
h = 10

bottomIntense = true
#= bottomIntense = false =#
if bottomIntense
    κ = @. κ0 + κ1*exp(-(ẑ + H0)/h)
else
    κ1 = 1e-4
    κ = κ1*ones(nẑ)
end

# timestepping
#= Δt = 86400 =#
Δt = 3600
nDaysSave = 1
adaptiveTimestep = false
α = 0.0

# Ekman layer
δ = sqrt(2*Pr*κ1/abs(f))

# start with far-field geostrophic flow
#= v0 = 0.01 =#        
#= v0 = -0.01 =#        
#= v0 = 0.1 =#        
v0 = -0.1        

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

log(ofile, @sprintf("nẑ = %d", nẑ))
log(ofile, @sprintf("H0 = %d m", H0))
log(ofile, @sprintf("Pr = %1.1f", Pr))
log(ofile, @sprintf("f  = %1.1e s-1", f))
log(ofile, @sprintf("N  = %1.1e s-1", N))
log(ofile, @sprintf("θ  = %1.1e rad", θ))
log(ofile, @sprintf("κ0 = %1.1e m2 s-1", κ0))
log(ofile, @sprintf("κ1 = %1.1e m2 s-1", κ1))
log(ofile, @sprintf("h  = %d m", h))
log(ofile, @sprintf("Δt = %.2f days", Δt/86400))
log(ofile, @sprintf("α = %.1f", α))
log(ofile, @sprintf("S = %.1f", S))

log(ofile, string("Canonical:              ", canonical))
log(ofile, string("Bottom intensification: ", bottomIntense))
log(ofile, string("Adaptive timestep:      ", adaptiveTimestep))

log(ofile, @sprintf("\nEkman layer thickness ~ %1.2f m", δ))
log(ofile, @sprintf("          ẑ[2] - ẑ[1] ~ %1.2f m", ẑ[2] - ẑ[1]))

log(ofile, @sprintf("\nSpin-down ~ %.1e days\n", abs(1/S/f)/86400))
close(ofile)
