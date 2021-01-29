# parameters

# canonical or transport-constrained case?
#= canonical = true =#
canonical = false

Ek = 1e-6 
S = 3.1e-4
Pr = 1e3
κ1 = 1e-3
v0 = -1

H = 1/sqrt(Ek) # z ∈ [0, H] ⟹ z̃ ∈ [0, H/δ = 1/sqrt(Ek)]
κ0 = 0
h = 0

# nondim arrest time
τ_A = 1/S

# nondim spindown time
τ_S = 1/sqrt(Ek)

# timestep
Δt = minimum([τ_S, τ_A])/100
tSave = τ_A

# number of grid points
nẑ = 2^8

# grid
ẑ = @. H*(1 - cos(pi*(0:nẑ-1)/(nẑ-1)))/2 # chebyshev (ẑ = 0 is bottom)

#= bottomIntense = true =#
bottomIntense = false
if bottomIntense
    κ = @. κ0 + κ1*exp(-ẑ/h)
else
    κ = κ1*ones(nẑ)
end

# timestepping
adaptiveTimestep = false
α = 0.5

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
log(ofile, @sprintf("H  = %1.5e", H))
log(ofile, @sprintf("Pr = %1.5e", Pr))
log(ofile, @sprintf("S  = %1.5e", S))
log(ofile, @sprintf("κ0 = %1.5e", κ0))
log(ofile, @sprintf("κ1 = %1.5e", κ1))
log(ofile, @sprintf("h  = %1.5e", h))
log(ofile, @sprintf("v0 = %1.5e", v0))
log(ofile, @sprintf("Δt = %1.5e", Δt))
log(ofile, @sprintf("α  = %1.5e", α))

log(ofile, string("\nCanonical:              ", canonical))
log(ofile, string("Bottom intensification: ", bottomIntense))
log(ofile, string("Adaptive timestep:      ", adaptiveTimestep))
log(ofile, string("Adaptive timestep:      ", adaptiveTimestep))

log(ofile, @sprintf("\nτ_A/τ_S  = %1.5e", τ_A/τ_S))

close(ofile)
