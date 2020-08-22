# parameters (as in MR91)
L = 0.4
H0 = 0.12
Pr = 1e0
f = 0.66
N = 3

# bowl
c = 2*L # chord length
R = c^2/(8*H0) + H0/2 # circle radius
H(r) = sqrt(R^2 - r^2) + H0 - R
Hr(r) = -r/sqrt(R^2 - r^2)

# number of grid points
nρ = 2^8 + 1
nσ = 2^8

# domain in terrain-following (ρ, σ) space
dρ = dx = L/(nρ + 1)
ρ = dρ:dρ:(L - dρ) # avoid spots where H = 0 or r = 0
σ = @. -(cos(pi*(0:nσ-1)/(nσ-1)) + 1)/2 # chebyshev 

# domain in physical (r, z) space (2D arrays)
r = repeat(ρ, 1, nσ)
z = repeat(σ', nρ, 1).*repeat(H.(ρ), 1, nσ)
dz = zeros(nρ, nσ)
dz[:, 1:end-1] = z[:, 2:end] - z[:, 1:end-1]
dz[:, end] = dz[:, end-1]

# diffusivity
κ0 = 1e-8
κ1 = 1e-6
h = H0/10
#= κ = κ1*ones(nρ, nσ) =#
κ = @. κ0 + κ1*exp(-(z + H(r))/h)
println("Ekman layer thickness ~ ", sqrt(2*Pr*κ1/abs(f)))
println("z[2] - z[1] ~ ", H0*(σ[2] - σ[1]))
