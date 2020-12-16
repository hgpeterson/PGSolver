# parameters (as in RC20)
H0 = 2e3
Pr = 1e0
f = -5.5e-5
N = 1e-3
θ = 2.5e-3

# canonical or transport-constrained case?
#= canonical = true =#
canonical = false

# number of grid points
nẑ = 2^8

# grid
#= z = @. -H0*(cos(pi*(0:nẑ-1)/(nẑ-1)) + 1)/2 # chebyshev =# 
#= ẑ = z/cos(θ) =#
ẑ = @. -H0*(cos(pi*(0:nẑ-1)/(nẑ-1)) + 1)/2 # chebyshev 

# diffusivity
κ0 = 6e-5
κ1 = 2e-3
h = 200
bottomIntense = true
if bottomIntense
    κ = @. κ0 + κ1*exp(-(ẑ + H0)/h)
else
    κ = κ1*ones(nẑ)
end

# timestepping
Δt = 86400
adaptiveTimestep = false
α = 0.0

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

log(ofile, string("Canonical:              ", canonical))
log(ofile, string("Bottom intensification: ", bottomIntense))
log(ofile, string("Adaptive timestep:      ", adaptiveTimestep))

log(ofile, @sprintf("\nEkman layer thickness ~ %1.2f m", sqrt(2*Pr*κ1/abs(f))))
log(ofile, @sprintf("          ẑ[2] - ẑ[1] ~ %1.2f m\n", ẑ[2] - ẑ[1]))
close(ofile)
