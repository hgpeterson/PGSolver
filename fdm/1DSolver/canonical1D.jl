################################################################################
# Utility functions 
################################################################################

"""
    u, w = rotate(û)

Rotate `û` into physical coordinate components `u` and `w`.
"""
function rotate(û)
    u = @. û*cosθ
    w = @. û*sinθ
    return u, w
end

################################################################################
# Solver functions 
################################################################################

"""
    diffMat, diffVec, bdyMat, explicitMat = getMatrices()   

Compute matrices for 1D equations.
"""
function getMatrices()
    nVars = 3
    nPts = nVars*nx*nz

    umap = reshape(1:nPts, nVars, nx, nz)    
    diffMat = Tuple{Int64,Int64,Float64}[]         # diffusion operator matrix 
    diffVec = zeros(nPts)                          # diffusion operator vector 
    bdyMat = Tuple{Int64,Int64,Float64}[]          # matrix for boundary conditions
    explicitMat = Tuple{Int64,Int64,Float64}[]     # matrix for explicit parts of right hand side

    # Main loop, insert stencil in matrices for each node point
    for i=1:nx
        for j=2:nz-1
            # dẑ stencil
            fd_ẑ = mkfdstencil(ẑẑ[i, j-1:j+1], ẑẑ[i, j], 1)
            κ_ẑ = sum(fd_ẑ.*κ[i, j-1:j+1])

            # dẑẑ stencil
            fd_ẑẑ = mkfdstencil(ẑẑ[i, j-1:j+1], ẑẑ[i, j], 2)

            # u diffusion term: dẑ(Pr*κ*dẑ(u)) = Pr*dẑ(κ)*dẑ(u) + Pr*κ*dẑẑ(u)
            row = umap[1, i, j]
            push!(diffMat, (row, umap[1, i, j-1], Pr*(κ_ẑ*fd_ẑ[1] + κ[i, j]*fd_ẑẑ[1])))
            push!(diffMat, (row, umap[1, i, j],   Pr*(κ_ẑ*fd_ẑ[2] + κ[i, j]*fd_ẑẑ[2])))
            push!(diffMat, (row, umap[1, i, j+1], Pr*(κ_ẑ*fd_ẑ[3] + κ[i, j]*fd_ẑẑ[3])))

            # v diffusion term: dẑ(Pr*κ*dẑ(v)) = Pr*dẑ(κ)*dẑ(v) + Pr*κ*dẑẑ(v)
            row = umap[2, i, j]
            push!(diffMat, (row, umap[2, i, j-1], Pr*(κ_ẑ*fd_ẑ[1] + κ[i, j]*fd_ẑẑ[1])))
            push!(diffMat, (row, umap[2, i, j],   Pr*(κ_ẑ*fd_ẑ[2] + κ[i, j]*fd_ẑẑ[2])))
            push!(diffMat, (row, umap[2, i, j+1], Pr*(κ_ẑ*fd_ẑ[3] + κ[i, j]*fd_ẑẑ[3])))

            # b diffusion term: dẑ(κ(N^2*cos(θ) + dẑ(b))) = dẑ(κ)*N^2*cos(θ) + dẑ(κ)*dẑ(b) + κ*dẑẑ(b)
            row = umap[3, i, j]
            push!(diffMat, (row, umap[3, i, j-1], (κ_ẑ*fd_ẑ[1] + κ[i, j]*fd_ẑẑ[1])))
            push!(diffMat, (row, umap[3, i, j],   (κ_ẑ*fd_ẑ[2] + κ[i, j]*fd_ẑẑ[2])))
            push!(diffMat, (row, umap[3, i, j+1], (κ_ẑ*fd_ẑ[3] + κ[i, j]*fd_ẑẑ[3])))
            diffVec[row] = κ_ẑ*N^2*cosθ[i, j]

            # 1st eqtn rhs: f*v*cos(θ) + b*sin(θ)
            row = umap[1, i, j]
            push!(explicitMat, (row, umap[2, i, j], f*cosθ[i, j]))
            push!(explicitMat, (row, umap[3, i, j], sinθ[i, j]))
            # 2nd eqtn rhs: -f*u*cos(θ)
            row = umap[2, i, j]
            push!(explicitMat, (row, umap[1, i, j], -f*cosθ[i, j]))
            # 3rd eqtn rhs: -u*N^2*sin(θ)
            row = umap[3, i, j]
            push!(explicitMat, (row, umap[1, i, j], -N^2*sinθ[i, j]))
        end

        # Boundary Conditions: Bottom
        # u = 0
        row = umap[1, i, 1] 
        push!(bdyMat, (row, row, 1.0))
        # v = 0
        row = umap[2, i, 1] 
        push!(bdyMat, (row, row, 1.0))
        # dẑ(b) = -N^2*cos(θ)
        row = umap[3, i, 1] 
        fd_ẑ = mkfdstencil(ẑẑ[i, 1:3], ẑẑ[i, 1], 1)
        push!(bdyMat, (row, umap[3, i, 1], fd_ẑ[1]))
        push!(bdyMat, (row, umap[3, i, 2], fd_ẑ[2]))
        push!(bdyMat, (row, umap[3, i, 3], fd_ẑ[3]))

        # Boundary Conditions: Top
        fd_ẑ = mkfdstencil(ẑẑ[i, nz-2:nz], ẑẑ[i, nz], 1)
        # dσ(u) = 0
        row = umap[1, i, nz] 
        push!(bdyMat, (row, umap[1, i, nz-2], fd_ẑ[1]))
        push!(bdyMat, (row, umap[1, i, nz-1], fd_ẑ[2]))
        push!(bdyMat, (row, umap[1, i, nz],   fd_ẑ[3]))
        # dσ(v) = 0
        row = umap[2, i, nz] 
        push!(bdyMat, (row, umap[2, i, nz-2], fd_ẑ[1]))
        push!(bdyMat, (row, umap[2, i, nz-1], fd_ẑ[2]))
        push!(bdyMat, (row, umap[2, i, nz],   fd_ẑ[3]))
        # dẑ(b) = 0
        row = umap[3, i, nz]
        push!(bdyMat, (row, umap[3, i, nz-2], fd_ẑ[1]))
        push!(bdyMat, (row, umap[3, i, nz-1], fd_ẑ[2]))
        push!(bdyMat, (row, umap[3, i, nz],   fd_ẑ[3]))
    end

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
    sol = evolveCanonical1D(nSteps)

Solve 1D equations over 2D ridge with time.
"""
function evolveCanonical1D(nSteps)
    # grid points
    nVars = 3
    nPts = nVars*nx*nz

    # timestep
    Δt = 3*3600
    nStepsPlot = 800
    nStepsSave = 80

    # for flattening for matrix mult
    umap = reshape(1:nPts, nVars, nx, nz)    
    bottomBdy = umap[:, :, 1][:]
    topBdy = umap[:, :, nz][:]

    # get matrices and vectors
    diffMat, diffVec, bdyMat, explicitMat = getMatrices()

    # left-hand side for evolution equation (save LU decomposition for speed)
    LHS = lu(getLHS(Δt, diffMat, bdyMat, bottomBdy, topBdy))

    # initial condition
    t = 0
    sol = zeros(nVars, nx, nz)
    #= # load data =#
    #= file = h5open("b.h5", "r") =#
    #= b = read(file, "b") =#
    #= t = read(file, "t") =#
    #= close(file) =#

    # plot initial state
    û = sol[1, :, :]
    v = sol[2, :, :]
    b = sol[3, :, :]
    iImg = 0
    plotCurrentState(t, û, v, b, iImg)

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
        RHS[umap[1, :, 1]]  .= 0 # u = 0 bot
        RHS[umap[1, :, nz]] .= 0 # u decay top
        RHS[umap[2, :, 1]]  .= 0 # v = 0 bot
        RHS[umap[2, :, nz]] .= 0 # v decay top
        RHS[umap[3, :, 1]]  = -N^2*cosθ[:, 1] # b flux bot
        RHS[umap[3, :, nz]] .= 0 # b flux top

        # solve
        solVec = LHS\RHS

        # log
        println(@sprintf("t = %.2f days (i = %d)", tDays, i))
        if i % nStepsPlot == 0
            # gather solution and rotate
            sol = reshape(solVec, 3, nx, nz)
            û = sol[1, :, :]
            v = sol[2, :, :]
            b = sol[3, :, :]

            # plot flow
            iImg += 1
            plotCurrentState(t, û, v, b, iImg)
        end
        if i % nStepsSave == 0
            # gather solution and rotate
            sol = reshape(solVec, 3, nx, nz)
            û = sol[1, :, :]
            v = sol[2, :, :]
            b = sol[3, :, :]

            # save data
            filename = @sprintf("sol%d.h5", tDays)
            println("saving ", filename)
            file = h5open(filename, "w")
            write(file, "û", û)
            write(file, "v", v)
            write(file, "b", b)
            write(file, "t", t)
            close(file)
        end
    end

    sol = reshape(solVec, 3, nx, nz)

    return sol
end

"""
    sol = steadyState()

Solve 1D equations over 2D ridge to steady state.
"""
function steadyState()
    # grid points
    nVars = 3
    nPts = nVars*nx*nz

    # for flattening for matrix mult
    umap = reshape(1:nPts, nVars, nx, nz)    
    bottomBdy = umap[:, :, 1][:]
    topBdy = umap[:, :, nz][:]

    # get matrices and vectors
    diffMat, diffVec, bdyMat, explicitMat = getMatrices()

    # LHS
    LHS = explicitMat + diffMat

    # boundaries
    LHS[bottomBdy, :] = bdyMat[bottomBdy, :]
    LHS[topBdy, :] = bdyMat[topBdy, :]

    # RHS
    RHS = -diffVec
    # boundaries
    RHS[umap[1, :, 1]]  .= 0 # u = 0 bot
    RHS[umap[1, :, nz]] .= 0 # u decay top
    RHS[umap[2, :, 1]]  .= 0 # v = 0 bot
    RHS[umap[2, :, nz]] .= 0 # v decay top
    RHS[umap[3, :, 1]]  .= -N^2*cosθ[:, 1] # b flux bot
    RHS[umap[3, :, nz]] .= 0    # b flux top

    # solve
    solVec = LHS\RHS

    # gather solution and rotate
    sol = reshape(solVec, 3, nx, nz)
    û = sol[1, :, :]
    v = sol[2, :, :]
    b = sol[3, :, :]

    # save data
    filename = "solSteady.h5"
    println("saving ", filename)
    file = h5open(filename, "w")
    write(file, "û", û)
    write(file, "v", v)
    write(file, "b", b)
    write(file, "t", Inf)
    close(file)

    # plot flow
    plotCurrentState(Inf, û, v, b, 999)

    return sol
end
