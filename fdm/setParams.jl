using Printf

# parameters (as in RC20)
L = 2e6
H0 = 2e3
Pr = 1e0
f = -5.5e-5
N = 1e-3

# topography
amp =  0.4*H0
H(x) = H0 - amp*sin(2*pi*x/L) # hill
Hx(x) = -2*pi/L*amp*cos(2*pi*x/L)

# number of grid points
nξ = 2^8 + 1 
nσ = 2^8

# domain in terrain-following (ξ, σ) space
dξ = dx = L/nξ
ξ = 0:dξ:(L - dξ)
σ = @. -(cos(pi*(0:nσ-1)/(nσ-1)) + 1)/2 # chebyshev 
ξξ = repeat(ξ, 1, nσ)
σσ = repeat(σ', nξ, 1)
dσ = zeros(nξ, nσ)
dσ[:, 1:end-1] = σσ[:, 2:end] - σσ[:, 1:end-1]
dσ[:, end] = dσ[:, end-1]

# domain in physical (x, z) space (2D arrays)
x = repeat(ξ, 1, nσ)
z = repeat(σ', nξ, 1).*repeat(H.(ξ), 1, nσ)

# arrays of sin(θ) and cos(θ) for 1D solutions
sinθ = @. -Hx(ξξ)/sqrt(1 + Hx(ξξ)^2)
cosθ = @. 1/sqrt(1 + Hx(ξξ)^2) 

# diffusivity
κ0 = 6e-5
κ1 = 2e-3
h = 200
#= κ = κ1*ones(nξ, nσ) =#
κ = @. κ0 + κ1*exp(-(z + H(x))/h)

# print properties
println("\nPGSolver with Parameters\n")

println(@sprintf("nξ = %d", nξ))
println(@sprintf("nσ = %d\n", nσ))

println(@sprintf("L  = %d km", L/1000))
println(@sprintf("H0 = %d m", H0))
println(@sprintf("Pr = %1.1f", Pr))
println(@sprintf("f  = %1.1e s-1", f))
println(@sprintf("N  = %1.1e s-1", N))
println(@sprintf("κ0 = %1.1e m2 s-1", κ0))
println(@sprintf("κ1 = %1.1e m2 s-1", κ1))
println(@sprintf("h  = %d m", h))

println(@sprintf("\nEkman layer thickness ~ %1.2f m", sqrt(2*Pr*κ1/abs(f))))
println(@sprintf("          z[2] - z[1] ~ %1.2f m", H0*(σ[2] - σ[1])))
