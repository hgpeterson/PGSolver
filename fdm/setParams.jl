# parameters (as in RC20)
L = 1e6
H0 = 1e3
Pr = 1e0
f = -5.5e-5
N = 1e-3

# topography
amp =  0.4*H0
#= amp =  0.04*H0 =#

#= H(x) = H0 # flat =#
#= Hx(x) = 0 =#

H(x) = H0 - amp*sin(2*pi*x/L) # hill
Hx(x) = -2*pi/L*amp*cos(2*pi*x/L)

# number of grid points
nξ = 2^8 + 1 
nσ = 2^8

# domain in terrain-following (ξ, σ) space
dξ = dx = L/nξ
ξ = 0:dξ:(L - dξ)
#= dσ = 1/(nσ - 1) =#
#= σ = -1:dσ:0 =#
σ = @. -(cos(pi*(0:nσ-1)/(nσ-1)) + 1)/2 # chebyshev 

# domain in physical (x, z) space (2D arrays)
x = repeat(ξ, 1, nσ)
z = repeat(σ', nξ, 1).*repeat(H.(ξ), 1, nσ)
dz = zeros(nξ, nσ)
dz[:, 1:end-1] = z[:, 2:end] - z[:, 1:end-1]
dz[:, end] = dz[:, end-1]

# diffusivity
κ0 = 6e-5
κ1 = 2e-3
h = 200
κ = κ1*ones(nξ, nσ)
#= κ = @. κ0 + κ1*exp(-(z + H(x))/h) =#
println("Ekman layer thickness ~ ", sqrt(2*Pr*κ1/abs(f)))
println("z[2] - z[1] ~ ", H0*(σ[2] - σ[1]))
