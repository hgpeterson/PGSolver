using SparseArrays, LinearAlgebra

include("myJuliaLib.jl")

"""
    diffMat, diffVec, bdyMat, explicitMat = getMatrices()   

Compute matrices for 1D equations.
"""
function getMatrices(κ, H, θ)
    nVars = 3
    nPts = nVars*nσ

    umap = reshape(1:nPts, nVars, nσ)    
    diffMat = Tuple{Int64,Int64,Float64}[]         # diffusion operator matrix 
    diffVec = zeros(nPts)                          # diffusion operator vector 
    bdyMat = Tuple{Int64,Int64,Float64}[]          # matrix for boundary conditions
    explicitMat = Tuple{Int64,Int64,Float64}[]     # matrix for explicit parts of right hand side

    # Main loop, insert stencil in matrices for each node point
    for j=2:nσ-1
        # dσ stencil
        fd_σ = mkfdstencil(σ[j-1:j+1], σ[j], 1)
        κ_σ = sum(fd_σ.*κ[j-1:j+1])

        # dσσ stencil
        fd_σσ = mkfdstencil(σ[j-1:j+1], σ[j], 2)

        # u diffusion term: dσ(Pr*κ*dσ(u)/H))/H = Pr*dσ(κ)*dσ(u)/H^2 + Pr*κ*dσσ(u)/H^2
        row = umap[1, j]
        push!(diffMat, (row, umap[1, j-1], Pr*(κ_σ*fd_σ[1] + κ[j]*fd_σσ[1])/H^2))
        push!(diffMat, (row, umap[1, j],   Pr*(κ_σ*fd_σ[2] + κ[j]*fd_σσ[2])/H^2))
        push!(diffMat, (row, umap[1, j+1], Pr*(κ_σ*fd_σ[3] + κ[j]*fd_σσ[3])/H^2))

        # v diffusion term: dσ(Pr*κ*dσ(v)/H))/H = Pr*dσ(κ)*dσ(v)/H^2 + Pr*κ*dσσ(v)/H^2
        row = umap[2, j]
        push!(diffMat, (row, umap[2, j-1], Pr*(κ_σ*fd_σ[1] + κ[j]*fd_σσ[1])/H^2))
        push!(diffMat, (row, umap[2, j],   Pr*(κ_σ*fd_σ[2] + κ[j]*fd_σσ[2])/H^2))
        push!(diffMat, (row, umap[2, j+1], Pr*(κ_σ*fd_σ[3] + κ[j]*fd_σσ[3])/H^2))

        # b diffusion term: dσ(κ(N^2 + dσ(b)/H))/H = N^2*cos(θ)*dσ(κ)/H + dσ(κ)*dσ(b)/H^2 + κ*dσσ(b)/H^2
        row = umap[3, j]
        push!(diffMat, (row, umap[3, j-1], (κ_σ*fd_σ[1] + κ[j]*fd_σσ[1])/H^2))
        push!(diffMat, (row, umap[3, j],   (κ_σ*fd_σ[2] + κ[j]*fd_σσ[2])/H^2))
        push!(diffMat, (row, umap[3, j+1], (κ_σ*fd_σ[3] + κ[j]*fd_σσ[3])/H^2))
        diffVec[row] = N^2*cos(θ)*κ_σ/H

        # 1st eqtn rhs: f*v*cos(θ) + b*sin(θ)
        row = umap[1, j]
        push!(explicitMat, (row, umap[2, j], f*cos(θ)))
        push!(explicitMat, (row, umap[3, j], sin(θ)))
        # 2nd eqtn rhs: -f*u*cos(θ)
        row = umap[2, j]
        push!(explicitMat, (row, umap[1, j], -f*cos(θ)))
        # 3rd eqtn rhs: -u*N^2*sin(θ)
        row = umap[3, j]
        push!(explicitMat, (row, umap[1, j], -N^2*sin(θ)))
    end

    # Boundary Conditions: Bottom
    # u = 0
    row = umap[1, 1] 
    push!(bdyMat, (row, row, 1.0))
    # v = 0
    row = umap[2, 1] 
    push!(bdyMat, (row, row, 1.0))
    # dσ(b)/H = -N^2*cos(θ)
    row = umap[3, 1] 
    fd_σ = mkfdstencil(σ[1:3], σ[1], 1)
    push!(bdyMat, (row, umap[3, 1], fd_σ[1]/H))
    push!(bdyMat, (row, umap[3, 2], fd_σ[2]/H))
    push!(bdyMat, (row, umap[3, 3], fd_σ[3]/H))

    # Boundary Conditions: Top
    fd_σ = mkfdstencil(σ[nσ-2:nσ], σ[nσ], 1)
    # dσ(u) = 0
    row = umap[1, nσ] 
    push!(bdyMat, (row, umap[1, nσ-2], fd_σ[1]/H))
    push!(bdyMat, (row, umap[1, nσ-1], fd_σ[2]/H))
    push!(bdyMat, (row, umap[1, nσ],   fd_σ[3]/H))
    # dσ(v) = 0
    row = umap[2, nσ] 
    push!(bdyMat, (row, umap[2, nσ-2], fd_σ[1]/H))
    push!(bdyMat, (row, umap[2, nσ-1], fd_σ[2]/H))
    push!(bdyMat, (row, umap[2, nσ],   fd_σ[3]/H))
    # dσ(b)/H = -N^2*cos(θ)
    row = umap[3, nσ]
    push!(bdyMat, (row, umap[3, nσ-2], fd_σ[1]/H))
    push!(bdyMat, (row, umap[3, nσ-1], fd_σ[2]/H))
    push!(bdyMat, (row, umap[3, nσ],   fd_σ[3]/H))

    # Create CSC sparse matrix from matrix elements
    diffMat = sparse((x->x[1]).(diffMat), (x->x[2]).(diffMat), (x->x[3]).(diffMat), nPts, nPts)
    bdyMat = sparse((x->x[1]).(bdyMat), (x->x[2]).(bdyMat), (x->x[3]).(bdyMat), nPts, nPts)
    explicitMat = sparse((x->x[1]).(explicitMat), (x->x[2]).(explicitMat), (x->x[3]).(explicitMat), nPts, nPts)

    return diffMat, diffVec, bdyMat, explicitMat
end

"""
    LHS = getLHS(Δt, diffMat, bdyMat, bottomBdy, topBdy)

Get implicit euler left hand side matrix.
"""
function getLHS(Δt, diffMat, bdyMat, bottomBdy, topBdy)
    # implicit euler
    LHS = I - diffMat*Δt 

    # no flux boundaries
    LHS[bottomBdy, :] = bdyMat[bottomBdy, :]
    LHS[topBdy, :] = bdyMat[topBdy, :]

    return LHS
end

"""
    u, v, w, b = rotate(û, v̂, b̂, θ)

Take the rotate variables and convert them to non-rotated fields.
"""
function rotate(û, v̂, b̂, θ)
    u = cos(θ)*û
    w = sin(θ)*û
    v = v̂
    b = b̂
    return u, v, w, b
end

"""
    sol = evolve1D(nSteps)

Solve the 1D equations at a point with time.
"""
function evolve1D(nSteps)
    # point on ridge
    iξ = 1
    Hcolumn = H(ξ[iξ])
    κcolumn = κ[iξ, :]
    θcolumn = acos(cosθ[iξ, 1])

    # grid points
    nVars = 3
    nPts = nVars*nσ

    # timestep
    Δt = 3*3600
    nStepsSave = 8000

    # for flattening for matrix mult
    umap = reshape(1:nPts, nVars, nσ)    
    bottomBdy = umap[:, 1]
    topBdy = umap[:, nσ]

    # get matrices and vectors
    diffMat, diffVec, bdyMat, explicitMat = getMatrices(κcolumn, Hcolumn, θcolumn)

    # left-hand side for evolution equation (save LU decomposition for speed)
    LHS = lu(getLHS(Δt, diffMat, bdyMat, bottomBdy, topBdy))

    # initial condition
    t = 0
    sol = zeros(nVars, nσ)

    # flatten for matrix mult
    solVec = reshape(sol, nPts, 1)

    # main loop
    for i=1:nSteps
        t += Δt
        tDays = t/86400

        # implicit euler diffusion
        diffRHS = solVec + diffVec*Δt

        # function to compute explicit RHS
        fExplicit(solVec, t) = explicitMat*solVec

        # explicit timestep for RHS
        explicitRHS = RK4(t, Δt, solVec, fExplicit)

        # sum the two
        RHS = diffRHS + explicitRHS

        # boundary conditions
        RHS[umap[1, 1]]  = 0 # u = 0 bot
        RHS[umap[1, nσ]] = 0 # u decay top
        RHS[umap[2, 1]]  = 0 # v = 0 bot
        RHS[umap[2, nσ]] = 0 # v decay top
        RHS[umap[3, 1]]  = -N^2*cos(θcolumn) # b flux bot
        #= RHS[umap[3, nσ]] = -N^2*cos(θcolumn) # b flux top =#
        RHS[umap[3, nσ]] = 0 # b flux top

        # solve
        solVec = LHS\RHS

        # log
        println(@sprintf("t = %.2f days (i = %d)", tDays, i))
        if i % nStepsSave == 0
            # gather solution and rotate
            sol = reshape(solVec, 3, nσ)
            û = sol[1, :]
            v̂ = sol[2, :]
            b̂ = sol[3, :]
            u, v, w, b = rotate(û, v̂, b̂, θcolumn)

            # save data
            file = h5open(@sprintf("sol%d.h5", tDays), "w")
            write(file, "u", u)
            write(file, "v", v)
            write(file, "w", w)
            write(file, "b", b)
            write(file, "t", t)
            close(file)
        end
    end

    sol = reshape(solVec, 3, nσ)

    return sol
end

"""
    sol = steadyState1D()

Solve the 1D equations at a point.
"""
function steadyState1D()
    # point on ridge
    iξ = 1
    Hcolumn = H(ξ[iξ])
    κcolumn = κ[iξ, :]
    θcolumn = acos(cosθ[iξ, 1])

    # grid points
    nVars = 3
    nPts = nVars*nσ

    # for flattening for matrix mult
    umap = reshape(1:nPts, nVars, nσ)    
    bottomBdy = umap[:, 1]
    topBdy = umap[:, nσ]

    # get matrices and vectors
    diffMat, diffVec, bdyMat, explicitMat = getMatrices(κcolumn, Hcolumn, θcolumn)

    # LHS
    LHS = explicitMat + diffMat
    # boundaries
    LHS[bottomBdy, :] = bdyMat[bottomBdy, :]
    LHS[topBdy, :] = bdyMat[topBdy, :]

    # RHS
    RHS = -diffVec
    # boundaries
    RHS[umap[1, 1]]  = 0 # u = 0 bot
    RHS[umap[1, nσ]] = 0 # u decay top
    RHS[umap[2, 1]]  = 0 # v = 0 bot
    RHS[umap[2, nσ]] = 0 # v decay top
    RHS[umap[3, 1]]  = -N^2*cos(θcolumn) # b flux bot
    RHS[umap[3, nσ]] = 0    # b flux top

    # solve
    solVec = LHS\RHS

    # gather solution and rotate
    sol = reshape(solVec, 3, nσ)
    û = sol[1, :]
    v̂ = sol[2, :]
    b̂ = sol[3, :]
    u, v, w, b = rotate(û, v̂, b̂, θcolumn)

    # save
    file = h5open("solSteady.h5", "w")
    write(file, "u", u)
    write(file, "v", v)
    write(file, "w", w)
    write(file, "b", b)
    write(file, "t", Inf)
    close(file)

    return sol
end
