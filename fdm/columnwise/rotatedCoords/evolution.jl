"""
    matrices = getEvolutionMatrices()

Compute the matrices needed for evolution equation integration.
"""
function getEvolutionMatrices()
    nPts = nx*nz

    umap = reshape(1:nPts, nx, nz)    
    diffMat = Tuple{Int64,Int64,Float64}[]         # diffusion operator matrix 
    diffVec = zeros(nPts)                          # diffusion operator vector 
    bdyFluxMat = Tuple{Int64,Int64,Float64}[]      # flux at boundary matrix

    # Main loop, insert stencil in matrices for each node point
    for i=1:nx
        # periodic in x
        iL = mod1(i-1, nx)
        iR = mod1(i+1, nx)

        # interior nodes only for operators
        for j=2:nz-1
            row = umap[i, j] 

            # dẑ stencil
            fd_ẑ = mkfdstencil(ẑ[i, j-1:j+1], ẑ[i, j], 1)
            κ_ẑ = sum(fd_ẑ.*κ[i, j-1:j+1])

            # dẑẑ stencil
            fd_ẑẑ = mkfdstencil(ẑ[i, j-1:j+1], ẑ[i, j], 2)

            # diffusion term: dẑ(κ(N^2*cos(θ) + dẑ(b))) = dẑ(κ)*N^2*cos(θ) + dẑ(κ)*dẑ(b) + κ*dẑẑ(b)
            push!(diffMat, (row, umap[i, j-1], (κ_ẑ*fd_ẑ[1] + κ[i, j]*fd_ẑẑ[1])))
            push!(diffMat, (row, umap[i, j],   (κ_ẑ*fd_ẑ[2] + κ[i, j]*fd_ẑẑ[2])))
            push!(diffMat, (row, umap[i, j+1], (κ_ẑ*fd_ẑ[3] + κ[i, j]*fd_ẑẑ[3])))
            diffVec[row] = κ_ẑ*N^2*cosθ[i, j]
        end

        # flux at boundaries: bottom
        row = umap[i, 1] 
        # dẑ stencil
        fd_ẑ = mkfdstencil(ẑ[i, 1:3], ẑ[i, 1], 1)
        # flux term: dẑ(b) = -N^2*cos(θ)
        push!(bdyFluxMat, (row, umap[i, 1], fd_ẑ[1]))
        push!(bdyFluxMat, (row, umap[i, 2], fd_ẑ[2]))
        push!(bdyFluxMat, (row, umap[i, 3], fd_ẑ[3]))

        # flux at boundaries: top
        row = umap[i, nz] 
        # dẑ stencil
        fd_ẑ = mkfdstencil(ẑ[i, nz-2:nz], ẑ[i, nz], 1)
        # flux term: dẑ(b) = 0
        push!(bdyFluxMat, (row, umap[i, nz-2], fd_ẑ[1]))
        push!(bdyFluxMat, (row, umap[i, nz-1], fd_ẑ[2]))
        push!(bdyFluxMat, (row, umap[i, nz],   fd_ẑ[3]))
    end

    # Create CSC sparse matrix from matrix elements
    diffMat = sparse((x->x[1]).(diffMat), (x->x[2]).(diffMat), (x->x[3]).(diffMat), nPts, nPts)
    bdyFluxMat = sparse((x->x[1]).(bdyFluxMat), (x->x[2]).(bdyFluxMat), (x->x[3]).(bdyFluxMat), nPts, nPts)

    return diffMat, diffVec, bdyFluxMat
end

"""
    evolutionLHS = getEvolutionLHS(Δt, diffMat, bdyFluxMat, bottomBdy, topBdy)

Generate the left-hand side matrix for the evolution problem of the form `I - diffmat*Δt`
and the not flux boundary condition applied to the boundaries
"""
function getEvolutionLHS(Δt, diffMat, bdyFluxMat, bottomBdy, topBdy)
    # implicit euler
    evolutionLHS = I - diffMat*Δt 

    # no flux boundaries
    evolutionLHS[bottomBdy, :] = bdyFluxMat[bottomBdy, :]
    evolutionLHS[topBdy, :] = bdyFluxMat[topBdy, :]

    return evolutionLHS
end

"""
    saveCheckpoint(b, chi, û, v, U, t)

Save .h5 checkpoint file for state `b` at time `t`.
"""
function saveCheckpoint(b, chi, û, v, U, t)
    tDays = t/86400
    savefile = @sprintf("checkpoint%d.h5", tDays)
    file = h5open(savefile, "w")
    write(file, "b", b)
    write(file, "chi", chi)
    write(file, "û", û)
    write(file, "v", v)
    write(file, "U", U)
    write(file, "t", t)
    write(file, "L", L)
    write(file, "H0", H0)
    write(file, "Pr", Pr)
    write(file, "f", f)
    write(file, "N", N)
    write(file, "symmetry", symmetry)
    write(file, "κ", κ)
    close(file)
    println(savefile)
end

"""
    b, chi, û, v, U, t, L, H0, Pr, f, N, symmetry, κ = loadCheckpoint(filename)

Load .h5 checkpoint file given by `filename`.
"""
function loadCheckpoint(filename)
    file = h5open(filename, "r")
    b = read(file, "b")
    chi = read(file, "chi")
    û = read(file, "û")
    v = read(file, "v")
    U = read(file, "U")
    t = read(file, "t")
    L = read(file, "L")
    H0 = read(file, "H0")
    Pr = read(file, "Pr")
    f = read(file, "f")
    N = read(file, "N")
    symmetry = read(file, "symmetry")
    κ = read(file, "κ")
    close(file)
    return b, chi, û, v, U, t, L, H0, Pr, f, N, symmetry, κ
end

"""
    b = evolve(tFinalDays)

Solve equation for `b` for `nSteps` time steps.
"""
function evolve(tFinalDays)
    # grid points
    nPts = nx*nz

    # timestep
    nSteps = Int64(tFinalDays*86400/Δt)
    nStepsInvert = 1
    nDaysPlot = 500
    nDaysSave = 10
    nStepsPlot = Int64(nDaysPlot*86400/Δt)
    nStepsSave = Int64(nDaysSave*86400/Δt)

    # for flattening for matrix mult
    umap = reshape(1:nPts, nx, nz)    
    bottomBdy = umap[:, 1]
    topBdy = umap[:, nz]

    # get matrices and vectors
    diffMat, diffVec, bdyFluxMat = getEvolutionMatrices()

    # left-hand side for evolution equation (save LU decomposition for speed)
    evolutionLHS = lu(getEvolutionLHS(Δt, diffMat, bdyFluxMat, bottomBdy, topBdy))

    # initial condition
    t = 0
    b = zeros(nx, nz)
    chi, û, v, U = invert(b)
    saveCheckpoint(b, chi, û, v, U, t)

    # flatten for matrix mult
    bVec = reshape(b, nPts, 1)
    ûVec = reshape(û, nPts, 1)
    sinθVec = reshape(sinθ, nPts, 1)
    
    # plot initial state of all zeros and no flow
    iImg = 0
    plotCurrentState(t, chi, û, v, b, iImg)

    # main loop
    for i=1:nSteps
        t += Δt
        tDays = t/86400

        # implicit euler diffusion
        diffRHS = bVec + diffVec*Δt

        # compute advection RHS
        advRHS = @. -ûVec*N^2*sinθVec

        # sum the two
        evolutionRHS = diffRHS + advRHS

        # boundary fluxes
        evolutionRHS[bottomBdy] .= -N^2*cosθ[:, 1] 
        evolutionRHS[topBdy] .= 0

        # solve
        bVec = evolutionLHS\evolutionRHS

        # log
        println(@sprintf("t = %.2f days (i = %d)", tDays, i))

        # reshape
        b = reshape(bVec, nx, nz)

        # invert buoyancy for flow
        chi, û, v, U = invert(b)
        ûVec = reshape(û, nPts, 1)

        # plot and/or save
        if i % nStepsPlot == 0
            iImg += 1
            plotCurrentState(t, chi, û, v, b, iImg)
        end
        if i % nStepsSave == 0
            saveCheckpoint(b, chi, û, v, U, t)
        end
    end

    b = reshape(bVec, nx, nz)

    return b
end
