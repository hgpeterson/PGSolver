################################################################################
# Functions used to compute the flow field given a buoyancy perturbation using
# finite differences, terrain-following coordinates, and taking advantage of 
# a 2D geometry.
################################################################################
using SparseArrays, LinearAlgebra, SpecialFunctions

include("myJuliaLib.jl")
include("terrainFollowing.jl")

"""
    LHS = getInversionLHS()

Setup left hand side of linear system for problem.
"""
function getInversionLHS()
    nPts = nξ*nσ
    iU = nPts + 1 # add equation for vertically integrated zonal flow

    umap = reshape(1:nPts, nξ, nσ)    
    A = Tuple{Int64,Int64,Float64}[]  

    # for finite difference on the top and bottom boundary
    fd_bot = mkfdstencil(σ[1:3], σ[1], 1)
    fd_top = mkfdstencil(σ[nσ-2:nσ], σ[nσ], 1)
    fd_top_σσ = mkfdstencil(σ[nσ-3:nσ], σ[nσ], 2)

    # Main loop, insert stencil in matrix for each node point
    for i=1:nξ
        # Boundary conditions: periodic in ξ
        iL = mod1(i-1, nξ)
        iR = mod1(i+1, nξ)

        # Lower boundary condition: chi = 0
        push!(A, (umap[i, 1], umap[i, 1], 1.0))

        # Upper boundary condition: chi = U
        push!(A, (umap[i, nσ], umap[i, nσ], 1.0))
        push!(A, (umap[i, nσ], iU, -1.0))

        # Interior nodes
        for j=2:nσ-1
            row = umap[i, j] 

            # dσσ stencil
            fd_σσ = mkfdstencil(σ[j-1:j+1], σ[j], 2)

            # eqtn: (f^2 + r^2)*dσσ(chi)/r/H^2 = -dξ(b) + dx(H)*σ*dσ(b)/H
            push!(A, (row, umap[i, j-1], (f^2 + r^2)*fd_σσ[1]/r/H(ξ[i])^2))
            push!(A, (row, umap[i, j],   (f^2 + r^2)*fd_σσ[2]/r/H(ξ[i])^2))
            push!(A, (row, umap[i, j+1], (f^2 + r^2)*fd_σσ[3]/r/H(ξ[i])^2))
        end
    end

    # just set U to 0
    row = iU
    push!(A, (row, iU, 1.0))

    # Create CSC sparse matrix from matrix elements
    A = sparse((x->x[1]).(A), (x->x[2]).(A), (x->x[3]).(A), nPts + 1, nPts + 1)

    return A
end

"""
    RHS = getInversionRHS(b)

Setup right hand side of linear system for problem.
"""
function getInversionRHS(b)
    nPts = nξ*nσ
    iU = nPts + 1 # add equation for vertically integrated zonal flow

    umap = reshape(1:nPts, nξ, nσ)    
    A = Tuple{Int64,Int64,Float64}[]  
    rhs = zeros(nPts + 1)                      

    # Main loop, insert stencil in matrix for each node point
    for i=1:nξ
        # Boundary conditions: periodic in ξ
        iL = mod1(i-1, nξ)
        iR = mod1(i+1, nξ)

        # Interior nodes
        for j=2:nσ-1
            row = umap[i, j] 

            # dσ stencil
            fd_σ = mkfdstencil(σ[j-1:j+1], σ[j], 1)
            
            # eqtn: (f^2 + r^2)*dσσ(chi)/r/H^2 = -dξ(b) + dx(H)*σ*dσ(b)/H
            rhs[row] = -(b[iR, j] - b[iL, j])/(2*dξ) + Hx(ξ[i])*σ[j]*sum(fd_σ.*b[i, j-1:j+1])/H(ξ[i])
        end
    end

    return rhs
end

"""
    sol = solveSystem(b)

Setup and solve linear system.
"""
function solveSystem(b)
    # get matrix and vector
    LHS = getInversionLHS()
    RHS = getInversionRHS(b)

    # solve 
    sol = LHS \ RHS
    return sol
end

"""
    chi, uξ, uη, uσ, U = postProcess(sol)

Take solution `sol` and extract reshaped `chi` and `U`. Compute `uξ`, `uη`, `uσ` 
from definition of chi.
"""
function postProcess(sol)
    # chi at σ = 0 is vertical integral of uξ
    U = sol[end] #FIXME think about this for TF coords

    # reshape rest of solution to get chi
    chi = reshape(sol[1:end-1], nξ, nσ)

    # compute uξ = dσ(chi)/H
    uξ = σDerivativeTF(chi)./H.(x)

    # compute uη = -f*uξ/r
    uη = -f*uξ/r

    # compute uσ = -dξ(chi)/H
    uσ = -ξDerivativeTF(chi)./H.(x)

    return chi, uξ, uη, uσ, U
end

"""
    chi, uξ, uη, uσ, U = invert(b, inversionLHS)

Wrapper function that inverts for flow given buoyancy perturbation `b`.
"""
function invert(b, inversionLHS)
    # compute RHS
    inversionRHS = getInversionRHS(b)

    # solve
    sol = inversionLHS\inversionRHS

    # compute flow from sol
    chi, uξ, uη, uσ, U = postProcess(sol)

    return chi, uξ, uη, uσ, U
end

"""
    b, chi, uξ, uη, uσ, U = pointwise1D(t, inversionLHS)

Apply the 1D solution to the Rayleigh drag problem pointwise over the domain.
See CF18 for details.
"""
function pointwise1D(t)
    # inverse boundary layer thickness
    q = @. sqrt(r*N^2*Hx(x)^2/(κ*(f^2 + r^2)))

    # time dependent analytical buoyancy solution (only works for constant κ)
    z1D = @. z + H(x)
    b = @. N^2/sqrt(1 + Hx(x)^2)/q*(exp(-q*z1D) - 0.5*(exp(-q*z1D)*erfc(q*sqrt(κ*t) - z1D/2/sqrt(κ*t)) + exp(q*z1D)*erfc(q*sqrt(κ*t) + z1D/2/sqrt(κ*t))))

    #TODO IN ORDER TO INVERT FOR FLOW, NEED ROTATED INVERSION EQTN
    return b
end
