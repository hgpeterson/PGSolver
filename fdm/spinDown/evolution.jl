"""
    A, diffVec = getRHS()   

If sol = (u, v, b, P_x) then we compute `A` and `diffVec` such that sol_t = A*sol + diffVec
and `A` also contains the proper boundary conditions.
"""
function getMatrices()
    nVars = 3
    nPts = nVars*nẑ + 1

    umap = reshape(1:(nPts-1), nVars, nẑ)    
    A = Tuple{Int64,Int64,Float64}[]  
    diffVec = zeros(nPts)

    # Main loop, insert stencil in matrices for each node point
    for j=2:nẑ-1
        # dẑ stencil
        fd_ẑ = mkfdstencil(ẑ[j-1:j+1], ẑ[j], 1)
        κ_ẑ = sum(fd_ẑ.*κ[j-1:j+1])

        # dẑẑ stencil
        fd_ẑẑ = mkfdstencil(ẑ[j-1:j+1], ẑ[j], 2)

        # 1st eqtn: u_t = f*v*cos(θ) - P_x*cos(θ) + b*sin(θ) + (nu*u_ẑ)_ẑ
        row = umap[1, j]
        # first term
        push!(A, (row, umap[2, j], f*cos(θ)))
        # second term
        push!(A, (row, nPts, -cos(θ)))
        # third term
        push!(A, (row, umap[3, j], sin(θ)))
        # fourth term: dẑ(Pr*κ*dẑ(u))) = Pr*dẑ(κ)*dẑ(u) + Pr*κ*dẑẑ(u)
        push!(A, (row, umap[1, j-1], Pr*(κ_ẑ*fd_ẑ[1] + κ[j]*fd_ẑẑ[1])))
        push!(A, (row, umap[1, j],   Pr*(κ_ẑ*fd_ẑ[2] + κ[j]*fd_ẑẑ[2])))
        push!(A, (row, umap[1, j+1], Pr*(κ_ẑ*fd_ẑ[3] + κ[j]*fd_ẑẑ[3])))

        # 2nd eqtn: v_t = -f*u*cos(θ) + (nu*v_ẑ)_ẑ
        row = umap[2, j]
        # first term:
        push!(A, (row, umap[1, j], -f*cos(θ)))
        # second term: dẑ(Pr*κ*dẑ(v))) = Pr*dẑ(κ)*dẑ(v) + Pr*κ*dẑẑ(v)
        push!(A, (row, umap[2, j-1], Pr*(κ_ẑ*fd_ẑ[1] + κ[j]*fd_ẑẑ[1])))
        push!(A, (row, umap[2, j],   Pr*(κ_ẑ*fd_ẑ[2] + κ[j]*fd_ẑẑ[2])))
        push!(A, (row, umap[2, j+1], Pr*(κ_ẑ*fd_ẑ[3] + κ[j]*fd_ẑẑ[3])))

        # 3rd eqtn: b_t = -u*N^2*sin(θ) + [κ*(N^2*cos(θ) + b_ẑ)]_ẑ
        row = umap[3, j]
        # first term
        push!(A, (row, umap[1, j], -N^2*sin(θ)))
        # second term: dẑ(κ(N^2*cos(θ) + dẑ(b))) = N^2*cos(θ)*dẑ(κ) + dẑ(κ)*dẑ(b) + κ*dẑẑ(b)
        push!(A, (row, umap[3, j-1], (κ_ẑ*fd_ẑ[1] + κ[j]*fd_ẑẑ[1])))
        push!(A, (row, umap[3, j],   (κ_ẑ*fd_ẑ[2] + κ[j]*fd_ẑẑ[2])))
        push!(A, (row, umap[3, j+1], (κ_ẑ*fd_ẑ[3] + κ[j]*fd_ẑẑ[3])))
        diffVec[row] = N^2*cos(θ)*κ_ẑ
    end

    # Boundary Conditions: Bottom
    # u = 0
    row = umap[1, 1] 
    push!(A, (row, row, 1.0))
    # v = 0
    row = umap[2, 1] 
    push!(A, (row, row, 1.0))
    # b_ẑ = -N^2*cos(θ)
    row = umap[3, 1] 
    fd_ẑ = mkfdstencil(ẑ[1:3], ẑ[1], 1)
    push!(A, (row, umap[3, 1], fd_ẑ[1]))
    push!(A, (row, umap[3, 2], fd_ẑ[2]))
    push!(A, (row, umap[3, 3], fd_ẑ[3]))

    # Boundary Conditions: Top
    fd_ẑ = mkfdstencil(ẑ[nẑ-2:nẑ], ẑ[nẑ], 1)
    # dẑ(u) = 0
    row = umap[1, nẑ] 
    push!(A, (row, umap[1, nẑ-2], fd_ẑ[1]))
    push!(A, (row, umap[1, nẑ-1], fd_ẑ[2]))
    push!(A, (row, umap[1, nẑ],   fd_ẑ[3]))
    # dẑ(v) = 0
    row = umap[2, nẑ] 
    push!(A, (row, umap[2, nẑ-2], fd_ẑ[1]))
    push!(A, (row, umap[2, nẑ-1], fd_ẑ[2]))
    push!(A, (row, umap[2, nẑ],   fd_ẑ[3]))
    # dẑ(b) = -N^2*cos(θ)
    row = umap[3, nẑ]
    push!(A, (row, umap[3, nẑ-2], fd_ẑ[1]))
    push!(A, (row, umap[3, nẑ-1], fd_ẑ[2]))
    push!(A, (row, umap[3, nẑ],   fd_ẑ[3]))

    # transport constraint
    row = nPts
    if canonical == true
        # canonical 1D: P_x = constant in time
        push!(A, (row, nPts, 1.0))
    else
        # transport-constrained 1D: P_x such that ∫u dẑ = 0
        for j=1:nẑ-1
            # trapeẑoidal rule: (u_j+1 + u_j)/2 * Δẑ_j
            push!(A, (row, umap[1, j],   (ẑ[j+1] - ẑ[j])/2))
            push!(A, (row, umap[1, j+1], (ẑ[j+1] - ẑ[j])/2))
        end
    end

    # Create CSC sparse matrix from matrix elements
    A = sparse((x->x[1]).(A), (x->x[2]).(A), (x->x[3]).(A), nPts, nPts)

    return A, diffVec
end

"""
    LHS, RHS = getLHSandRHS(Δt, A, α, bottomBdy, topBdy)

Get left- and right-hand-side matrices for time stepping where
    (x^(n+1) - x^n)/Δt = α*A*x^n + (1 - α)*A*x^(n+1) + y
So that 
    LHS = I/Δt - (1 - α)*A,
    RHS = I/Δt + α*A.
"""
function getLHSandRHS(Δt, A, α, bottomBdy, topBdy)
    nVars = 3
    nPts = nVars*nẑ + 1

    # (x^(n+1) - x^n)/Δt = α*A*x^n + (1 - α)*A*x^(n+1)
    LHS = I/Δt - (1 - α)*A 
    RHS = I/Δt + α*A 

    # boundary conditions and integral constraint
    LHS[bottomBdy, :] = A[bottomBdy, :]
    LHS[topBdy, :] = A[topBdy, :]
    LHS[nPts, :] = A[nPts, :]

    return LHS, RHS
end

"""
    sol = evolve(tFinalDays)

Solve 1D equations for `tFinalDays` days.
"""
function evolve(tFinalDays)
    # grid points
    nVars = 3
    nPts = nVars*nẑ + 1

    # timestep
    nSteps = Int64(tFinalDays*86400/Δt)
    nStepsSave = Int64(nDaysSave*86400/Δt)

    # for flattening for matrix mult
    umap = reshape(1:(nPts-1), nVars, nẑ)    
    bottomBdy = umap[:, 1]
    topBdy = umap[:, nẑ]

    # get matrices and vectors
    A, diffVec = getMatrices()

    # left- and right-hand-side matrices for evolution
    LHS, RHS = getLHSandRHS(Δt, A, α, bottomBdy, topBdy)
    LHS = lu(LHS)

    # initial condition
    t = 0
    sol = zeros(nPts)
    # start with far-field geostrophic flow
    sol[umap[2, :]] .= v0
    sol[nPts] = f*v0
    # save initial condition
    û = sol[umap[1, :]]
    v = sol[umap[2, :]]
    b = sol[umap[3, :]]
    Px = sol[nPts]
    saveCheckpointSpinDown(û, v, b, Px, t)

    # main loop
    for i=1:nSteps
        t += Δt
        tDays = t/86400

        # right-hand-side as a vector
        RHSVec = RHS*sol + diffVec

        # boundary conditions
        RHSVec[umap[1, 1]]  = 0           # u = 0 bot
        RHSVec[umap[1, nẑ]] = 0           # u decay top
        RHSVec[umap[2, 1]]  = 0           # v = 0 bot
        RHSVec[umap[2, nẑ]] = 0           # v decay top
        RHSVec[umap[3, 1]]  = -N^2*cos(θ) # b flux bot
        #= RHSVec[umap[3, nẑ]] = -N^2*cos(θ) # b flux top =#
        RHSVec[umap[3, nẑ]] = 0           # b flux top
        if canonical
            RHSVec[nPts] = f*v0           # f*v0 = P_x
        else
            RHSVec[nPts] = 0              # integral constraint
        end

        # solve
        sol = LHS\RHSVec

        # log
        if i % nStepsSave == 0
            println(@sprintf("t = %.2f days (i = %d)", tDays, i))

            # gather solution
            û = sol[umap[1, :]]
            v = sol[umap[2, :]]
            b = sol[umap[3, :]]
            Px = sol[nPts]

            # save data
            saveCheckpointSpinDown(û, v, b, Px, t)
        end
    end

    return sol
end

#= """ =#
#=     sol = steadyState1D() =#

#= Solve the 1D equations at a point. =#
#= """ =#
#= function steadyState1D() =#
#=     # point on ridge =#
#=     iξ = 1 =#
#=     Hcolumn = H(ξ[iξ]) =#
#=     κcolumn = κ[iξ, :] =#
#=     θcolumn = acos(cosθ[iξ, 1]) =#

#=     # grid points =#
#=     nVars = 3 =#
#=     nPts = nVars*nẑ =#

#=     # for flattening for matrix mult =#
#=     umap = reshape(1:nPts, nVars, nẑ) =#    
#=     bottomBdy = umap[:, 1] =#
#=     topBdy = umap[:, nẑ] =#

#=     # get matrices and vectors =#
#=     diffMat, diffVec, bdyMat, explicitMat = getMatrices(κcolumn, Hcolumn, θcolumn) =#

#=     # LHS =#
#=     LHS = explicitMat + diffMat =#
#=     # boundaries =#
#=     LHS[bottomBdy, :] = bdyMat[bottomBdy, :] =#
#=     LHS[topBdy, :] = bdyMat[topBdy, :] =#

#=     # RHS =#
#=     RHS = -diffVec =#
#=     # boundaries =#
#=     RHS[umap[1, 1]]  = 0 # u = 0 bot =#
#=     RHS[umap[1, nẑ]] = 0 # u decay top =#
#=     RHS[umap[2, 1]]  = 0 # v = 0 bot =#
#=     RHS[umap[2, nẑ]] = 0 # v decay top =#
#=     RHS[umap[3, 1]]  = -N^2*cos(θcolumn) # b flux bot =#
#=     RHS[umap[3, nẑ]] = 0    # b flux top =#

#=     # solve =#
#=     solVec = LHS\RHS =#

#=     # gather solution and rotate =#
#=     sol = reshape(solVec, 3, nẑ) =#
#=     û = sol[1, :] =#
#=     v̂ = sol[2, :] =#
#=     b̂ = sol[3, :] =#
#=     u, v, w, b = rotate(û, v̂, b̂, θcolumn) =#

#=     # save =#
#=     file = h5open("solSteady.h5", "w") =#
#=     write(file, "u", u) =#
#=     write(file, "v", v) =#
#=     write(file, "w", w) =#
#=     write(file, "b", b) =#
#=     write(file, "t", Inf) =#
#=     close(file) =#

#=     return sol =#
#= end =#
