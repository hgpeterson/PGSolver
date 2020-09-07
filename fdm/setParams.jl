# parameters (as in RC20)
L = 4e6
H0 = 4e3
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
#= dz = zeros(nξ, nσ) =#
#= dz[:, 1:end-1] = z[:, 2:end] - z[:, 1:end-1] =#
#= dz[:, end] = dz[:, end-1] =#

# diffusivity
κ0 = 6e-5
κ1 = 2e-3
h = 200
#= κ = κ1*ones(nξ, nσ) =#
κ = @. κ0 + κ1*exp(-(z + H(x))/h)
println("Ekman layer thickness ~ ", sqrt(2*Pr*κ1/abs(f)))
println("z[2] - z[1] ~ ", H0*(σ[2] - σ[1]))
