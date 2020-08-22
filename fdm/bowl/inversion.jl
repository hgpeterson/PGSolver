################################################################################
# Functions used to compute the flow field given a buoyancy perturbation using
# finite differences, terrain-following coordinates, and taking advantage of 
# a 2D geometry.
################################################################################
using SparseArrays, PyPlot, LinearAlgebra

close("all")
plt.style.use("~/presentation_plots.mplstyle")

include("setParams.jl")
include("myJuliaLib.jl")
include("plottingLib.jl")

"""
    LHS = getInversionLHS()

Setup left hand side of linear system for problem.
"""
function getInversionLHS()
    nPts = nρ*nσ
    iU = nPts + 1 # add equation for vertically integrated zonal flow

    umap = reshape(1:nPts, nρ, nσ)    
    A = Tuple{Int64,Int64,Float64}[]  

    # for finite difference on the top and bottom boundary
    fd_bot = mkfdstencil(σ[1:3], σ[1], 1)
    fd_top = mkfdstencil(σ[nσ-2:nσ], σ[nσ], 1)
    fd_top_σσ = mkfdstencil(σ[nσ-3:nσ], σ[nσ], 2)

    
    # Left and Right boundary conditions: no flow at the center and edge
    for j=1:nσ
        push!(A, (umap[1, j], umap[1, j], 1.0))
        push!(A, (umap[nρ, j], umap[nρ, j], 1.0))
    end

    # the rest of the interior
    for i=2:nρ-1
        # Lower boundary conditions 
        # b.c. 1: dσ(chi) = 0
        push!(A, (umap[i, 1], umap[i, 1], fd_bot[1]))
        push!(A, (umap[i, 1], umap[i, 2], fd_bot[2]))
        push!(A, (umap[i, 1], umap[i, 3], fd_bot[3]))
        # b.c. 2: chi = 0 
        push!(A, (umap[i, 2], umap[i, 1], 1.0))

        # Upper boundary conditions
        # b.c. 1: dσσ(chi) = 0 (or -stress at top)
        push!(A, (umap[i, nσ], umap[i, nσ-3], fd_top_σσ[1]))
        push!(A, (umap[i, nσ], umap[i, nσ-2], fd_top_σσ[2]))
        push!(A, (umap[i, nσ], umap[i, nσ-1], fd_top_σσ[3]))
        push!(A, (umap[i, nσ], umap[i, nσ],   fd_top_σσ[4]))
        #= rhs[umap[1, i, nσ]] = -0.1 # some wind stress at the top =#
        # b.c. 2: chi = bottom stress 
        push!(A, (umap[i, nσ-1], umap[i, nσ], 1.0))
        push!(A, (umap[i, nσ-1], iU, -1.0))

        # Interior nodes
        for j=3:nσ-2
            row = umap[i, j] 

            # dσ stencil
            fd_σ = mkfdstencil(σ[j-1:j+1], σ[j], 1)
            κ_σ = sum(fd_σ.*κ[i, j-1:j+1])

            # dσσ stencil
            fd_σσ = mkfdstencil(σ[j-1:j+1], σ[j], 2)
            κ_σσ = sum(fd_σσ.*κ[i, j-1:j+1])

            # dσσσ stencil
            fd_σσσ = mkfdstencil(σ[j-2:j+2], σ[j], 3)

            # dσσσσ stencil
            fd_σσσσ = mkfdstencil(σ[j-2:j+2], σ[j], 4)
            
            # eqtn: dσσ(nu*dσσ(chi))/H^4 + f^2*chi/nu + f*U = dρ(b) - dx(H)*σ*dσ(b)/H
            # term 1 (product rule)
            push!(A, (row, umap[i, j-1], Pr*κ_σσ*fd_σσ[1]/H(ρ[i])^4))
            push!(A, (row, umap[i, j],   Pr*κ_σσ*fd_σσ[2]/H(ρ[i])^4))
            push!(A, (row, umap[i, j+1], Pr*κ_σσ*fd_σσ[3]/H(ρ[i])^4))

            push!(A, (row, umap[i, j-2], 2*Pr*κ_σ*fd_σσσ[1]/H(ρ[i])^4))
            push!(A, (row, umap[i, j-1], 2*Pr*κ_σ*fd_σσσ[2]/H(ρ[i])^4))
            push!(A, (row, umap[i, j],   2*Pr*κ_σ*fd_σσσ[3]/H(ρ[i])^4))
            push!(A, (row, umap[i, j+1], 2*Pr*κ_σ*fd_σσσ[4]/H(ρ[i])^4))
            push!(A, (row, umap[i, j+2], 2*Pr*κ_σ*fd_σσσ[5]/H(ρ[i])^4))

            push!(A, (row, umap[i, j-2], Pr*κ[i, j]*fd_σσσσ[1]/H(ρ[i])^4))
            push!(A, (row, umap[i, j-1], Pr*κ[i, j]*fd_σσσσ[2]/H(ρ[i])^4))
            push!(A, (row, umap[i, j],   Pr*κ[i, j]*fd_σσσσ[3]/H(ρ[i])^4))
            push!(A, (row, umap[i, j+1], Pr*κ[i, j]*fd_σσσσ[4]/H(ρ[i])^4))
            push!(A, (row, umap[i, j+2], Pr*κ[i, j]*fd_σσσσ[5]/H(ρ[i])^4))
            # term 2
            push!(A, (row, umap[i, j], f^2/(Pr*κ[i, j])))
            # term 3
            push!(A, (row, iU, f))
            # term 4 (rhs)
        end
    end

    # zonal mean / integral equation for bottom stress
    #   < int_-1^0 (f^2*chi/nu + U) dσ + dσ(nu*dσσ(chi))/H^3 > = 0
    row = iU
    push!(A, (row, iU, 1.0))
    #= rhs[row] = κ1/Hx(ρ[1]) =#
    
    #= # term 1 =#
    #= for i=1:nρ =#
    #=     for j=1:nσ =#
    #=         push!(A, (row, umap[i, j], f^2/(Pr*κ)*dρ/L*dσ)) =#
    #=         push!(A, (row, iU,         f*dρ/L*dσ)) =#
    #=     end =#
    #= end =#
    #= # term 2 =#
    #= fd_σσσ_top = mkfdstencil(σ[nσ-4:nσ], σ[nσ], 3) =#
    #= for i=1:nρ =#
    #=     push!(A, (row, umap[i, nσ-4], Pr*κ*fd_σσσ_top[1]*dρ/L/H(ρ[i])^3)) =#
    #=     push!(A, (row, umap[i, nσ-3], Pr*κ*fd_σσσ_top[2]*dρ/L/H(ρ[i])^3)) =#
    #=     push!(A, (row, umap[i, nσ-2], Pr*κ*fd_σσσ_top[3]*dρ/L/H(ρ[i])^3)) =#
    #=     push!(A, (row, umap[i, nσ-1], Pr*κ*fd_σσσ_top[4]*dρ/L/H(ρ[i])^3)) =#
    #=     push!(A, (row, umap[i, nσ],   Pr*κ*fd_σσσ_top[5]*dρ/L/H(ρ[i])^3)) =#
    #= end =#

    # Create CSC sparse matrix from matrix elements
    A = sparse((x->x[1]).(A), (x->x[2]).(A), (x->x[3]).(A), nPts + 1, nPts + 1)

    #= println("rank(A)  = ", rank(A)) =#
    #= println("nPts + 1 = ", nPts + 1) =#

    return A
end

"""
    RHS = getInversionRHS(b)

Setup right hand side of linear system for problem.
"""
function getInversionRHS(b)
    nPts = nρ*nσ
    iU = nPts + 1 # add equation for vertically integrated zonal flow

    umap = reshape(1:nPts, nρ, nσ)    
    A = Tuple{Int64,Int64,Float64}[]  
    rhs = zeros(nPts + 1)                      

    # Main loop, insert stencil in matrix for each node point
    for i=2:nρ-1
        iL = i-1
        iR = i+1

        # Interior nodes
        for j=3:nσ-2
            row = umap[i, j] 

            # dσ stencil
            fd_σ = mkfdstencil(σ[j-1:j+1], σ[j], 1)
            
            # eqtn: dσσ(nu*dσσ(chi))/H^4 + f^2*chi/nu + f*U = dρ(b) - dρ(H)*σ*dσ(b)/H
            rhs[row] = (b[iR, j] - b[iL, j])/(2*dρ) - Hr(ρ[i])*σ[j]*sum(fd_σ.*b[i, j-1:j+1])/H(ρ[i])
        end
    end

    # zonal mean / integral equation for bottom stress
    #   < int_-1^0 (f^2*chi/nu + U) dσ + dσ(nu*dσσ(chi))/H^3 > = 0
    
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
    chi, u, v, w, U = postProcess(sol)

Take solution `sol` and extract reshaped `chi` and `U`. Compute `u`, `v`, `w` 
from definition of chi.
"""
function postProcess(sol)
    # chi at σ = 0 is vertical integral of zonal flow
    U = sol[end]

    # reshape rest of solution to get chi
    chi = reshape(sol[1:end-1], nρ, nσ)

    # compute u = dz(chi)
    u = zDerivativeTF(chi)

    # compute v = int_-H^0 f*chi/nu dz = int_-1^0 f*chi/nu dσ*H
    v = zeros(nρ, nσ)
    for i=1:nρ
        v[i, :] = cumtrapz(f*chi[i, :]./(Pr*κ[i, :]) .- U, σ)*H(ρ[i])
    end

    # compute w = -dr(r*chi)/r 
    w = -rDerivativeTF(r.*chi)./r

    return chi, u, v, w, U
end

"""
    chi, u, v, w, U = invert(b, inversionLHS)

Wrapper function that inverts for flow given buoyancy perturbation `b`.
"""
function invert(b, inversionLHS)
    # compute RHS
    inversionRHS = getInversionRHS(b)

    # solve
    sol = inversionLHS\inversionRHS

    # compute flow from sol
    chi, u, v, w, U = postProcess(sol)

    return chi, u, v, w, U
end

"""
    chiEkman = getChiEkman(br)

Compute Ekman layer solution to problem given forcing dr(b).
"""
function getChiEkman(br)
    # Ekman layer thickness
    δ = sqrt(2*Pr*κ1/abs(f)) # using κ at the bottom

    # interior solution: thermal wind balance
    chi_I = br
    chi_I_bot = repeat(chi_I[:, 1], 1, nσ)
    chi_I_top = repeat(chi_I[:, nσ], 1, nσ)

    # bottom Ekman layer correction
    chi_B_bot = @. -exp(-(z + H(r))/δ)*chi_I_bot*(cos((z + H(r))/δ) + sin((z + H(r))/δ))

    # top Ekman layer correction
    chi_B_top = @. -exp(z/δ)*chi_I_top*cos(z/δ)

    # full solution (use full κ with assumption that its variation is larger than δ)
    chiEkman = @. κ/f^2*(chi_I + chi_B_bot + chi_B_top)

    return chiEkman
end

#= b = @. h*N^2*exp(-(z + H(r))/h) =#
#= br = ρDerivativeTF(b) =#
#= #1= br = @. -Hr(r)/h*b =1# =#

#= inversionLHS = getInversionLHS() =#
#= chi, u, v, w, U = invert(b, inversionLHS) =#
#= chiEkman = getChiEkman(br) =#
#= plotCurrentState(0, chi, chiEkman, u, v, w, b, 0) =#
