################################################################################
# Functions used to compute the flow field given a buoyancy perturbation using
# finite differences, terrain-following coordinates, and taking advantage of 
# a 2D geometry.
################################################################################

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
            
            # eqtn: dσσ(nu*dσσ(chi))/H^4 + f^2*(chi - U)/nu = dξ(b) - dx(H)*σ*dσ(b)/H
            # term 1 (product rule)
            push!(A, (row, umap[i, j-1], Pr*κ_σσ*fd_σσ[1]/H(ξ[i])^4))
            push!(A, (row, umap[i, j],   Pr*κ_σσ*fd_σσ[2]/H(ξ[i])^4))
            push!(A, (row, umap[i, j+1], Pr*κ_σσ*fd_σσ[3]/H(ξ[i])^4))

            push!(A, (row, umap[i, j-2], 2*Pr*κ_σ*fd_σσσ[1]/H(ξ[i])^4))
            push!(A, (row, umap[i, j-1], 2*Pr*κ_σ*fd_σσσ[2]/H(ξ[i])^4))
            push!(A, (row, umap[i, j],   2*Pr*κ_σ*fd_σσσ[3]/H(ξ[i])^4))
            push!(A, (row, umap[i, j+1], 2*Pr*κ_σ*fd_σσσ[4]/H(ξ[i])^4))
            push!(A, (row, umap[i, j+2], 2*Pr*κ_σ*fd_σσσ[5]/H(ξ[i])^4))

            push!(A, (row, umap[i, j-2], Pr*κ[i, j]*fd_σσσσ[1]/H(ξ[i])^4))
            push!(A, (row, umap[i, j-1], Pr*κ[i, j]*fd_σσσσ[2]/H(ξ[i])^4))
            push!(A, (row, umap[i, j],   Pr*κ[i, j]*fd_σσσσ[3]/H(ξ[i])^4))
            push!(A, (row, umap[i, j+1], Pr*κ[i, j]*fd_σσσσ[4]/H(ξ[i])^4))
            push!(A, (row, umap[i, j+2], Pr*κ[i, j]*fd_σσσσ[5]/H(ξ[i])^4))
            # term 2
            push!(A, (row, umap[i, j], f^2/(Pr*κ[i, j])))
            # term 3
            push!(A, (row, iU, f^2/(Pr*κ[i, j])))
            # term 4 (rhs)
        end
    end

    # zonal mean / integral equation for bottom stress
    #   < int_-1^0 (f^2*chi/nu + U) dσ + dσ(nu*dσσ(chi))/H^3 > = 0
    row = iU
    push!(A, (row, iU, 1.0))
    #= rhs[row] = κ1/Hx(ξ[1]) =#
    
    #= # term 1 =#
    #= for i=1:nξ =#
    #=     for j=1:nσ =#
    #=         push!(A, (row, umap[i, j], f^2/(Pr*κ)*dξ/L*dσ)) =#
    #=         push!(A, (row, iU,         f*dξ/L*dσ)) =#
    #=     end =#
    #= end =#
    #= # term 2 =#
    #= fd_σσσ_top = mkfdstencil(σ[nσ-4:nσ], σ[nσ], 3) =#
    #= for i=1:nξ =#
    #=     push!(A, (row, umap[i, nσ-4], Pr*κ*fd_σσσ_top[1]*dξ/L/H(ξ[i])^3)) =#
    #=     push!(A, (row, umap[i, nσ-3], Pr*κ*fd_σσσ_top[2]*dξ/L/H(ξ[i])^3)) =#
    #=     push!(A, (row, umap[i, nσ-2], Pr*κ*fd_σσσ_top[3]*dξ/L/H(ξ[i])^3)) =#
    #=     push!(A, (row, umap[i, nσ-1], Pr*κ*fd_σσσ_top[4]*dξ/L/H(ξ[i])^3)) =#
    #=     push!(A, (row, umap[i, nσ],   Pr*κ*fd_σσσ_top[5]*dξ/L/H(ξ[i])^3)) =#
    #= end =#

    # Create CSC sparse matrix from matrix elements
    A = sparse((x->x[1]).(A), (x->x[2]).(A), (x->x[3]).(A), nPts + 1, nPts + 1)

    return A
end

"""
    RHS = getInversionRHS(b)

Setup right hand side of linear system for problem.
"""
function getInversionRHS(b; ξVariation=true)
    nPts = nξ*nσ
    umap = reshape(1:nPts, nξ, nσ)    

    # eqtn: dσσ(nu*dσσ(chi))/H^4 + f^2*chi/nu + f*U = dξ(b) - dx(H)*σ*dσ(b)/H
    if ξVariation
        rhs = xDerivativeTF(b)
    else
        rhs = -Hx.(ξξ).*σσ.*σDerivativeTF(b)./H.(ξξ)
    end
    rhs[umap[:, [1, 2, nσ-2, nσ]]] .= 0 # boundary conditions require zeros on the rhs

    # return as vector
    rhsVec = zeros(nPts + 1)                      
    rhsVec[umap[:, :]] = reshape(rhs, nPts, 1)

    return rhsVec
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

    # compute uη = int_-1^0 f*chi/nu dσ*H
    uη = zeros(nξ, nσ)
    for i=1:nξ
        uη[i, :] = cumtrapz(f*chi[i, :]./(Pr*κ[i, :]) .- U, σ)*H(ξ[i])
    end

    # compute uσ = -dξ(chi)/H
    uσ = -ξDerivativeTF(chi)./H.(x)

    return chi, uξ, uη, uσ, U
end

"""
    chi, uξ, uη, uσ, U = invert(b, inversionLHS)

Wrapper function that inverts for flow given buoyancy perturbation `b`.
"""
function invert(b, inversionLHS; ξVariation=true)
    # compute RHS
    inversionRHS = getInversionRHS(b; ξVariation=ξVariation)

    # solve
    sol = inversionLHS\inversionRHS

    # compute flow from sol
    chi, uξ, uη, uσ, U = postProcess(sol)

    return chi, uξ, uη, uσ, U
end

"""
    chiEkman = getChiEkman(b)

Compute Ekman layer solution to problem given buoyancy perturbation b.
"""
function getChiEkman(b)
    # compute x derivative of b
    bx = xDerivativeTF(b)

    # Ekman layer thickness
    δ = sqrt(2*Pr*κ1/abs(f)) # using κ at the bottom

    # interior solution: thermal wind balance
    chi_I = bx
    chi_I_bot = repeat(chi_I[:, 1], 1, nσ)
    chi_I_top = repeat(chi_I[:, nσ], 1, nσ)

    # bottom Ekman layer correction
    chi_B_bot = @. -exp(-(z + H(x))/δ)*chi_I_bot*(cos((z + H(x))/δ) + sin((z + H(x))/δ))

    # top Ekman layer correction
    chi_B_top = @. -exp(z/δ)*chi_I_top*cos(z/δ)

    # full solution (use full κ with assumption that its variation is larger than δ)
    chiEkman = @. κ/f^2*(chi_I + chi_B_bot + chi_B_top)

    return chiEkman
end

#= """ =#
#=     LHS = getInversionLHS1DAdjusted() =#

#= Setup left hand side of linear system for adjusted 1D problem. =#
#= """ =#
#= function getInversionLHS1DAdjusted() =#
#=     nPts = nξ*nσ =#

#=     umap = reshape(1:nPts, nξ, nσ) =#    
#=     A = Tuple{Int64,Int64,Float64}[] =#  

#=     # for finite difference on the top and bottom boundary =#
#=     fd_bot = mkfdstencil(σ[1:3], σ[1], 1) =#
#=     fd_top = mkfdstencil(σ[nσ-2:nσ], σ[nσ], 1) =#
#=     fd_top_σσ = mkfdstencil(σ[nσ-3:nσ], σ[nσ], 2) =#

#=     # Main loop, insert stencil in matrix for each node point =#
#=     for i=1:nξ =#
#=         # Lower boundary conditions =# 
#=         # b.c. 1: dσ(chi) = 0 =#
#=         push!(A, (umap[i, 1], umap[i, 1], fd_bot[1])) =#
#=         push!(A, (umap[i, 1], umap[i, 2], fd_bot[2])) =#
#=         push!(A, (umap[i, 1], umap[i, 3], fd_bot[3])) =#
#=         # b.c. 2: chi = 0 =# 
#=         push!(A, (umap[i, 2], umap[i, 1], 1.0)) =#

#=         # Upper boundary conditions =#
#=         # b.c. 1: dσ(chi) = 0 =#
#=         push!(A, (umap[i, nσ], umap[i, nσ-2], fd_top[1])) =#
#=         push!(A, (umap[i, nσ], umap[i, nσ-1], fd_top[2])) =#
#=         push!(A, (umap[i, nσ], umap[i, nσ],   fd_top[3])) =#
#=         #1= # b.c. 1: dσσ(chi) = 0 (or -stress at top) =1# =#
#=         #1= push!(A, (umap[i, nσ], umap[i, nσ-3], fd_top_σσ[1])) =1# =#
#=         #1= push!(A, (umap[i, nσ], umap[i, nσ-2], fd_top_σσ[2])) =1# =#
#=         #1= push!(A, (umap[i, nσ], umap[i, nσ-1], fd_top_σσ[3])) =1# =#
#=         #1= push!(A, (umap[i, nσ], umap[i, nσ],   fd_top_σσ[4])) =1# =#
#=         # b.c. 2: chi = U =#
#=         push!(A, (umap[i, nσ-1], umap[i, nσ], 1.0)) =#

#=         # Interior nodes =#
#=         for j=3:nσ-2 =#
#=             row = umap[i, j] =# 

#=             # dσ stencil =#
#=             fd_σ = mkfdstencil(σ[j-1:j+1], σ[j], 1) =#
#=             κ_σ = sum(fd_σ.*κ[i, j-1:j+1]) =#

#=             # dσσ stencil =#
#=             fd_σσ = mkfdstencil(σ[j-1:j+1], σ[j], 2) =#
#=             κ_σσ = sum(fd_σσ.*κ[i, j-1:j+1]) =#

#=             # dσσσ stencil =#
#=             fd_σσσ = mkfdstencil(σ[j-2:j+2], σ[j], 3) =#

#=             # dσσσσ stencil =#
#=             fd_σσσσ = mkfdstencil(σ[j-2:j+2], σ[j], 4) =#
            
#=             # eqtn: dσσ(nu*dσσ(chi))/H^4 + f^2*cos(θ)^2(chi - U)/nu = -dσ(b)*sin(θ)/H =#
#=             # term 1 (product rule) =#
#=             push!(A, (row, umap[i, j-1], Pr*κ_σσ*fd_σσ[1]/H(ξ[i])^4)) =#
#=             push!(A, (row, umap[i, j],   Pr*κ_σσ*fd_σσ[2]/H(ξ[i])^4)) =#
#=             push!(A, (row, umap[i, j+1], Pr*κ_σσ*fd_σσ[3]/H(ξ[i])^4)) =#

#=             push!(A, (row, umap[i, j-2], 2*Pr*κ_σ*fd_σσσ[1]/H(ξ[i])^4)) =#
#=             push!(A, (row, umap[i, j-1], 2*Pr*κ_σ*fd_σσσ[2]/H(ξ[i])^4)) =#
#=             push!(A, (row, umap[i, j],   2*Pr*κ_σ*fd_σσσ[3]/H(ξ[i])^4)) =#
#=             push!(A, (row, umap[i, j+1], 2*Pr*κ_σ*fd_σσσ[4]/H(ξ[i])^4)) =#
#=             push!(A, (row, umap[i, j+2], 2*Pr*κ_σ*fd_σσσ[5]/H(ξ[i])^4)) =#

#=             push!(A, (row, umap[i, j-2], Pr*κ[i, j]*fd_σσσσ[1]/H(ξ[i])^4)) =#
#=             push!(A, (row, umap[i, j-1], Pr*κ[i, j]*fd_σσσσ[2]/H(ξ[i])^4)) =#
#=             push!(A, (row, umap[i, j],   Pr*κ[i, j]*fd_σσσσ[3]/H(ξ[i])^4)) =#
#=             push!(A, (row, umap[i, j+1], Pr*κ[i, j]*fd_σσσσ[4]/H(ξ[i])^4)) =#
#=             push!(A, (row, umap[i, j+2], Pr*κ[i, j]*fd_σσσσ[5]/H(ξ[i])^4)) =#
#=             # term 2 =# 
#=             push!(A, (row, umap[i, j], f^2*cosθ[i, j]^2/(Pr*κ[i, j]))) =#
#=         end =#
#=     end =#

#=     # Create CSC sparse matrix from matrix elements =#
#=     A = sparse((x->x[1]).(A), (x->x[2]).(A), (x->x[3]).(A), nPts, nPts) =#

#=     return A =#
#= end =#

#= """ =#
#=     RHS = getInversionRHS1DAdjusted(b) =#

#= Setup right hand side of linear system for the adjusted 1D problem. =#
#= """ =#
#= function getInversionRHS1DAdjusted(b) =#
#=     U = zeros(size(z)) # set vertically integrated cross-slope flow to zero =#

#=     nPts = nξ*nσ =#
#=     umap = reshape(1:nPts, nξ, nσ) =#    

#=     # -dẑ(b)*sin(θ) = -(-sin(θ)*dx(b) + cos(θ)*dz(b))*sin(θ) =#
#=     rhs = -(-sinθ.*xDerivativeTF(b) .+ cosθ.*zDerivativeTF(b)).*sinθ =#

#=     # reshape =#
#=     rhs = reshape(rhs, nPts, 1) =#

#=     # chi = U at top =#
#=     rhs[umap[:, nσ - 1]] = U[:, 1] =#

#=     return rhs =#
#= end =#

#= """ =#
#=     chi, u, v, w, U = postProcess1DAdjusted(sol) =#

#= Take solution `sol` for adjusted 1D problem and extract reshaped `chi`. Compute `u`, `v`, `w` =#
#= from definition of chi and equations. =#
#= """ =#
#= function postProcess1DAdjusted(sol) =#
#=     U = zeros(size(z)) # set vertically integrated cross-slope flow to zero =#

#=     # reshape =# 
#=     chi = reshape(sol, nξ, nσ) =#

#=     # compute û = dẑ(chi) =#
#=     û = -sinθ.*xDerivativeTF(chi) .+ cosθ.*zDerivativeTF(chi) =#

#=     # compute u =#
#=     u = cosθ.*û =#

#=     # compute w =#
#=     w = sinθ.*û =#

#=     # compute v = int_-1^0 f*chi/nu dσ*H =#
#=     v = zeros(nξ, nσ) =#
#=     for i=1:nξ =#
#=         v[i, :] = cumtrapz(f*cosθ[i, :].*(chi[i, :] .- U[i, :])./(Pr*κ[i, :]), σ)*H(ξ[i]) =#
#=     end =#
    
#=     return chi, u, v, w, U =#
#= end =#

#= """ =#
#=     chi, uξ, uη, uσ, U = invert1DAdjusted(b, inversionLHS1DAdjusted) =#

#= Wrapper function that inverts for flow given buoyancy perturbation `b`. =#
#= """ =#
#= function invert1DAdjusted(b, inversionLHS1DAdjusted) =#
#=     # compute RHS =#
#=     inversionRHS1DAdjusted = getInversionRHS1DAdjusted(b) =#

#=     # solve =#
#=     sol = inversionLHS1DAdjusted\inversionRHS1DAdjusted =#

#=     # compute flow from sol =#
#=     chi, u, v, w, U = postProcess1DAdjusted(sol) =#

#=     #1= # fixed v? =1# =#
#=     #1= Px = zeros((nξ, nσ)) =1# =#
#=     #1= bottomTurb = zDerivativeTF(Pr*κ.*zDerivativeTF(zDerivativeTF(chi))) =1# =#
#=     #1= for i=1:nξ =1# =#
#=     #1=     Px[i, :] .= b[i, 1]*sinθ[i, 1] + bottomTurb[i, 1] =1# =#
#=     #1= end =1# =#
#=     #1= v = @. -1/(f*cosθ)*(-Px + b*sinθ + bottomTurb) =1# =#

#=     # to TF =#
#=     uξ, uη, uσ = transformToTF(u, v, w) =#

#=     return chi, uξ, uη, uσ, U =#
#= end =#
