################################################################################
# Set parameters for PG model
################################################################################

L = 1e6
H0 = 1e3
amp = 0.4*H0

# hill
H(x) = H0 - amp*sin(2*pi*x/L - pi/2)
Hx(x) = -2*pi/L*amp*cos(2*pi*x/L - pi/2)

# flat
#= H(x) = H0 =# 
#= Hx(x) = 0 =# 

nu = 1e0
#= nu(x, z) = 1e0 =#
#= nu(x, z) = 1e2 + 1e4*exp(-(z + H(x))/200) =#
Pr = 1
f = 1e-4
N = 1e-3

println("Ekman layer thickness: ", sqrt(2*nu/f))
#= println("Ekman layer thickness: ", sqrt(2*nu(0, -H0)/f)) =#
