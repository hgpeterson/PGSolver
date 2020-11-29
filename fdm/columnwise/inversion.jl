################################################################################
# Functions used to compute the flow field given a buoyancy perturbation using
# finite differences, terrain-following coordinates, and taking advantage of 
# a 2D geometry.
################################################################################

"""
    LHS = getInversionLHS()

Setup left hand side of linear system for problem.
"""
function getInversionLHS(κ, H)
    iU = nσ + 1
    A = Tuple{Int64,Int64,Float64}[]  

    # for finite difference on the top and bottom boundary
    fd_bot = mkfdstencil(σ[1:3], σ[1], 1)
    fd_top = mkfdstencil(σ[nσ-2:nσ], σ[nσ], 1)
    fd_top_σσ = mkfdstencil(σ[nσ-3:nσ], σ[nσ], 2)

    # Main loop, insert stencil in matrix for each node point
    # Lower boundary conditions 
    # b.c. 1: dσ(chi) = 0
    push!(A, (1, 1, fd_bot[1]))
    push!(A, (1, 2, fd_bot[2]))
    push!(A, (1, 3, fd_bot[3]))
    # b.c. 2: chi = 0 
    push!(A, (2, 1, 1.0))

    # Upper boundary conditions
    # b.c. 1: dσσ(chi) = 0 
    push!(A, (nσ, nσ-3, fd_top_σσ[1]))
    push!(A, (nσ, nσ-2, fd_top_σσ[2]))
    push!(A, (nσ, nσ-1, fd_top_σσ[3]))
    push!(A, (nσ, nσ,   fd_top_σσ[4]))
    # b.c. 2: chi - U = 0
    push!(A, (nσ-1, nσ,  1.0))
    push!(A, (nσ-1, iU, -1.0))

    # Interior nodes
    for j=3:nσ-2
        row = j

        # dσ stencil
        fd_σ = mkfdstencil(σ[j-1:j+1], σ[j], 1)
        κ_σ = sum(fd_σ.*κ[j-1:j+1])

        # dσσ stencil
        fd_σσ = mkfdstencil(σ[j-1:j+1], σ[j], 2)
        κ_σσ = sum(fd_σσ.*κ[j-1:j+1])

        # dσσσ stencil
        fd_σσσ = mkfdstencil(σ[j-2:j+2], σ[j], 3)

        # dσσσσ stencil
        fd_σσσσ = mkfdstencil(σ[j-2:j+2], σ[j], 4)
        
        # eqtn: dσσ(nu*dσσ(chi))/H^4 + f^2*(chi - U)/nu = dξ(b) - dx(H)*σ*dσ(b)/H
        # term 1 (product rule)
        push!(A, (row, j-1, Pr*κ_σσ*fd_σσ[1]/H^4))
        push!(A, (row, j,   Pr*κ_σσ*fd_σσ[2]/H^4))
        push!(A, (row, j+1, Pr*κ_σσ*fd_σσ[3]/H^4))

        push!(A, (row, j-2, 2*Pr*κ_σ*fd_σσσ[1]/H^4))
        push!(A, (row, j-1, 2*Pr*κ_σ*fd_σσσ[2]/H^4))
        push!(A, (row, j,   2*Pr*κ_σ*fd_σσσ[3]/H^4))
        push!(A, (row, j+1, 2*Pr*κ_σ*fd_σσσ[4]/H^4))
        push!(A, (row, j+2, 2*Pr*κ_σ*fd_σσσ[5]/H^4))

        push!(A, (row, j-2, Pr*κ[j]*fd_σσσσ[1]/H^4))
        push!(A, (row, j-1, Pr*κ[j]*fd_σσσσ[2]/H^4))
        push!(A, (row, j,   Pr*κ[j]*fd_σσσσ[3]/H^4))
        push!(A, (row, j+1, Pr*κ[j]*fd_σσσσ[4]/H^4))
        push!(A, (row, j+2, Pr*κ[j]*fd_σσσσ[5]/H^4))
        # term 2
        push!(A, (row, j,   f^2/(Pr*κ[j])))
        push!(A, (row, iU, -f^2/(Pr*κ[j])))
    end

    # if dξ(b) ~ 0 then 
    #   (1) U = 0
    #       for fixed 1D solution
    #   (2) dσ(nu*dσσ(chi))/H^3 = Hx*b at σ = -1
    #       for canonical 1D solution
    row = iU
    if symmetry
        push!(A, (row, row, 1.0))
    else
        # dσ stencil
        fd_σ = mkfdstencil(σ[1:3], σ[1], 1)
        κ_σ = sum(fd_σ.*κ[1:3])

        # dσσ stencil
        fd_σσ = mkfdstencil(σ[1:4], σ[1], 2)

        # dσσσ stencil
        fd_σσσ = mkfdstencil(σ[1:5], σ[1], 3)

        # product rule
        push!(A, (row, 1, Pr*κ_σ*fd_σσ[1]/H^3))
        push!(A, (row, 2, Pr*κ_σ*fd_σσ[2]/H^3))
        push!(A, (row, 3, Pr*κ_σ*fd_σσ[3]/H^3))
        push!(A, (row, 4, Pr*κ_σ*fd_σσ[4]/H^3))

        push!(A, (row, 1, Pr*κ[1]*fd_σσσ[1]/H^3))
        push!(A, (row, 2, Pr*κ[1]*fd_σσσ[2]/H^3))
        push!(A, (row, 3, Pr*κ[1]*fd_σσσ[3]/H^3))
        push!(A, (row, 4, Pr*κ[1]*fd_σσσ[4]/H^3))
        push!(A, (row, 5, Pr*κ[1]*fd_σσσ[5]/H^3))
    end

    # Create CSC sparse matrix from matrix elements
    A = sparse((x->x[1]).(A), (x->x[2]).(A), (x->x[3]).(A), nσ+1, nσ+1)

    #= println(rank(A)) =#
    #= println(nσ+1) =#

    return A
end

"""
    RHS = getInversionRHS(b)

Setup right hand side of linear system for problem.
"""
function getInversionRHS(b)
    # last row is for U
    rhs = zeros(nξ, nσ+1)
    iU = nσ + 1

    # eqtn: dσσ(nu*dσσ(chi))/H^4 + f^2*(chi - U)/nu = dξ(b) - dx(H)*σ*dσ(b)/H
    if ξVariation
        rhs[:, 1:nσ] = xDerivativeTF(b)
    else
        rhs[:, 1:nσ] = -Hx.(ξξ).*σσ.*σDerivativeTF(b)./H.(ξξ)
    end
    rhs[:, [1, 2, nσ-1, nσ]] .= 0 # boundary conditions require zeros on the rhs

    # if dξ(b) ~ 0 then 
    #   (1) U = 0
    #       for fixed 1D solution
    #   (2) dσ(nu*dσσ(chi))/H^3 = Hx*b at σ = -1
    #       for canonical 1D solution
    if symmetry
        #= rhs[:, iU] .= 1.0 =#
    else
        rhs[:, iU] = Hx.(ξ).*b[:, 1]
    end

    return rhs
end

"""
    chi, uξ, uη, uσ, U = postProcess(sol)

Take solution `sol` and extract reshaped `chi` and `U`. Compute `uξ`, `uη`, `uσ` 
from definition of chi.
"""
function postProcess(sol)
    iU = nσ + 1

    # chi at σ = 0 is vertical integral of uξ
    U = sol[:, iU] 

    # rest of solution is chi
    chi = sol[:, 1:nσ]

    #= println(U[194]) =#
    #= println(chi[194, nσ]) =#

    # compute uξ = dσ(chi)/H
    uξ = σDerivativeTF(chi)./H.(x)

    # compute uη = int_-1^0 f*chi/nu dσ*H
    uη = zeros(nξ, nσ)
    for i=1:nξ
        uη[i, :] = cumtrapz(f*(chi[i, :] .- U[i])./(Pr*κ[i, :]), σ)*H(ξ[i])
    end

    # compute uσ = -dξ(chi)/H
    if ξVariation
        uσ = -ξDerivativeTF(chi)./H.(x)
    else
        uσ = zeros(nξ, nσ)
    end

    return chi, uξ, uη, uσ, U
end

"""
    chi, uξ, uη, uσ, U = invert(b)

Wrapper function that inverts for flow given buoyancy perturbation `b`.
"""
function invert(b)
    # compute RHS
    inversionRHS = getInversionRHS(b)

    # solve
    sol = zeros(nξ, nσ+1)
    for i=1:nξ
        sol[i, :] = inversionLHSs[i]\inversionRHS[i, :]
    end

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
