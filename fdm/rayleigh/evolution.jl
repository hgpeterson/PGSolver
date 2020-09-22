using SparseArrays, PyPlot, LinearAlgebra, Printf, HDF5

close("all")
plt.style.use("~/presentation_plots.mplstyle")
pygui(false)

include("setParams.jl")
include("myJuliaLib.jl")
include("plottingLib.jl")
include("inversion.jl")

"""
    matrices = getEvolutionMatrices()

Compute the matrices needed for evolution equation integration.
"""
function getEvolutionMatrices()
    nPts = nξ*nσ

    umap = reshape(1:nPts, nξ, nσ)    
    diffMat = Tuple{Int64,Int64,Float64}[]         # diffusion operator matrix 
    diffVec = zeros(nPts)                          # diffusion operator vector 
    bdyFluxMat = Tuple{Int64,Int64,Float64}[]      # flux at boundary matrix
    ξDerivativeMat = Tuple{Int64,Int64,Float64}[]  # advection operator matrix (ξ)
    σDerivativeMat = Tuple{Int64,Int64,Float64}[]  # advection operator matrix (σ)

    # Main loop, insert stencil in matrices for each node point
    for i=1:nξ
        # periodic in ξ
        iL = mod1(i-1, nξ)
        iR = mod1(i+1, nξ)

        # interior nodes only for operators
        for j=2:nσ-1
            row = umap[i, j] 

            # dσ stencil
            fd_σ = mkfdstencil(σ[j-1:j+1], σ[j], 1)
            κ_σ = sum(fd_σ.*κ[i, j-1:j+1])

            # dσσ stencil
            fd_σσ = mkfdstencil(σ[j-1:j+1], σ[j], 2)

            # diffusion term: dσ(κ(N^2 + dσ(b)/H))/H = N^2*dσ(κ)/H + dσ(κ)*dσ(b)/H^2 + κ*dσσ(b)/H^2
            push!(diffMat, (row, umap[i, j-1], (κ_σ*fd_σ[1] + κ[i, j]*fd_σσ[1])/H(ξ[i])^2))
            push!(diffMat, (row, umap[i, j],   (κ_σ*fd_σ[2] + κ[i, j]*fd_σσ[2])/H(ξ[i])^2))
            push!(diffMat, (row, umap[i, j+1], (κ_σ*fd_σ[3] + κ[i, j]*fd_σσ[3])/H(ξ[i])^2))
            diffVec[row] = N^2*κ_σ/H(ξ[i])

            # ξ advection term: dξ()
            push!(ξDerivativeMat, (row, umap[iR, j],  1.0/(2*dξ)))
            push!(ξDerivativeMat, (row, umap[iL, j], -1.0/(2*dξ)))

            # σ advection term: dσ()
            push!(σDerivativeMat, (row, umap[i, j-1], fd_σ[1]))
            push!(σDerivativeMat, (row, umap[i, j],   fd_σ[2]))
            push!(σDerivativeMat, (row, umap[i, j+1], fd_σ[3]))
        end

        # flux at boundaries: bottom
        row = umap[i, 1] 
        # dσ stencil
        fd_σ = mkfdstencil(σ[1:3], σ[1], 1)
        # flux term: dσ(b)/H = -N^2
        push!(bdyFluxMat, (row, umap[i, 1], fd_σ[1]/H(ξ[i])))
        push!(bdyFluxMat, (row, umap[i, 2], fd_σ[2]/H(ξ[i])))
        push!(bdyFluxMat, (row, umap[i, 3], fd_σ[3]/H(ξ[i])))

        # flux at boundaries: top
        row = umap[i, nσ] 
        # dσ stencil
        fd_σ = mkfdstencil(σ[nσ-2:nσ], σ[nσ], 1)
        # flux term: dσ(b)/H = -N^2
        push!(bdyFluxMat, (row, umap[i, nσ-2], fd_σ[1]/H(ξ[i])))
        push!(bdyFluxMat, (row, umap[i, nσ-1], fd_σ[2]/H(ξ[i])))
        push!(bdyFluxMat, (row, umap[i, nσ],   fd_σ[3]/H(ξ[i])))
    end

    # Create CSC sparse matrix from matrix elements
    diffMat = sparse((x->x[1]).(diffMat), (x->x[2]).(diffMat), (x->x[3]).(diffMat), nPts, nPts)
    bdyFluxMat = sparse((x->x[1]).(bdyFluxMat), (x->x[2]).(bdyFluxMat), (x->x[3]).(bdyFluxMat), nPts, nPts)
    ξDerivativeMat = sparse((x->x[1]).(ξDerivativeMat), (x->x[2]).(ξDerivativeMat), (x->x[3]).(ξDerivativeMat), nPts, nPts)
    σDerivativeMat = sparse((x->x[1]).(σDerivativeMat), (x->x[2]).(σDerivativeMat), (x->x[3]).(σDerivativeMat), nPts, nPts)

    return diffMat, diffVec, bdyFluxMat, ξDerivativeMat, σDerivativeMat
end

"""
    dy = explicitRHS(Δt, y, f)

Compute `dy`, the change in `y` given a timestep of `Δt` and that `dt(y) = f(y)`.
"""
function explicitRHS(Δt, y, f)
    #= # euler =#
    #= return Δt*f(y) =#
    #= # midpoint =#
	#= f1 = f(y) =#
    #= f2 = f(y + Δt*f1/2) =#
	#= dy = Δt*f2 =#
    #= return dy =#
    # RK4
	f1 = f(y)
    f2 = f(y + Δt*f1/2)
    f3 = f(y + Δt*f2/2)
    f4 = f(y + Δt*f3)
	dy = Δt/6*(f1 + 2*f2 + 2*f3 + f4)
    return dy
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
    b = evolveFullNL(nSteps)

Solve full nonlinear equation for `b` for `nSteps` time steps.
"""
function evolveFullNL(nSteps)
    # grid points
    nPts = nξ*nσ

    # timestep
    #= Δt = 43200 =#
    Δt = 10*86400
    nStepsInvert = 1
    nStepsPlot = 10
    nStepsProfile = 100
    adaptiveTimestep = false

    # for flattening for matrix mult
    umap = reshape(1:nPts, nξ, nσ)    
    bottomBdy = umap[:, 1]
    topBdy = umap[:, nσ]

    # get matrices and vectors
    diffMat, diffVec, bdyFluxMat, ξDerivativeMat, σDerivativeMat = getEvolutionMatrices()

    # left-hand side for evolution equation (save LU decomposition for speed)
    evolutionLHS = lu(getEvolutionLHS(Δt, diffMat, bdyFluxMat, bottomBdy, topBdy))

    # left-hand side for inversion equations
    inversionLHS = lu(getInversionLHS())

    # vector of H values
    HVec = reshape(H.(x), nPts, 1)

    # initial condition
    t = 0
    b = zeros(nξ, nσ)
    #= t = 800*86400 =#
    #= b, chi, uξ, uη, uσ, U = pointwise1D(800*86400, inversionLHS) =#
    #= # load data =#
    #= file = h5open("b.h5", "r") =#
    #= b = read(file, "b") =#
    #= t = read(file, "t") =#

    # invert initial condition
    chi, uξ, uη, uσ, U = invert(b, inversionLHS)
    uξVec = reshape(uξ, nPts, 1)
    uσVec = reshape(uσ, nPts, 1)
    
    # plot initial state of all zeros and no flow
    iImg = 0
    plotCurrentState(t, chi, uξ, uη, uσ, b, iImg)

    # initialize plot of profiles
    iξ = 1
    ax = profilePlotInit()

    # flatten for matrix mult
    bVec = reshape(b, nPts, 1)
    uξVec = reshape(uξ, nPts, 1)
    uσVec = reshape(uσ, nPts, 1)

    # main loop
    for i=1:nSteps
        t += Δt
        tDays = t/86400

        # implicit euler diffusion
        diffRHS = bVec + diffVec*Δt

        # function to compute advection RHS
        # (note the parentheses here to allow for sparse matrices to work first)
        #= fAdvRHS(bVec) = -(uξVec.*(ξDerivativeMat*bVec) + uσVec.*(σDerivativeMat*bVec) + uσVec.*HVec*N^2) =#
        #= println(maximum(abs.(fAdvRHS(bVec)))) =#
        
        # explicit timestep for advection
        #= advRHS = explicitRHS(Δt, bVec, fAdvRHS) =#
        #= println(maximum(abs.(advRHS))) =#
        #= bVec += advRHS =#

        # sum the two
        #= evolutionRHS = diffRHS + advRHS =#
        evolutionRHS = diffRHS 

        # no flux boundaries
        evolutionRHS[bottomBdy] .= -N^2
        evolutionRHS[topBdy] .= -N^2

        # solve
        bVec = evolutionLHS\evolutionRHS

        # log
        println(@sprintf("t = %.2f days (i = %d)", tDays, i))
        if i % nStepsInvert == 0
            # reshape
            b = reshape(bVec, nξ, nσ)

            # invert buoyancy for flow
            chi, uξ, uη, uσ, U = invert(b, inversionLHS)
            uξVec = reshape(uξ, nPts, 1)
            uσVec = reshape(uσ, nPts, 1)
            uξCFL = minimum(abs.(dξ./uξ))
            uσCFL = minimum(abs.(dσ./uσ))
            println(@sprintf("CFL uξ: %.2f days", uξCFL/86400))
            println(@sprintf("CFL uσ: %.2f days", uσCFL/86400))
            #= if 0.01*minimum([uξCFL, uσCFL]) < Δt && adaptiveTimestep =#
            #=     # need to have smaller step size by CFL =#
            #=     Δt = 0.01*minimum([uξCFL, uσCFL]) =#
            #=     println(@sprintf("Decreasing timestep to %.2f days", Δt/86400)) =#
            #=     evolutionLHS = lu(getEvolutionLHS(Δt, diffMat, bdyFluxMat, bottomBdy, topBdy)) =#
            #= elseif 0.01*minimum([uξCFL, uσCFL]) > 10*Δt && Δt < 10*86400 && adaptiveTimestep =#
            #=     # could have much larger step size by CFL =#
            #=     Δt = 0.01*minimum([uξCFL, uσCFL]) =#
            #=     println(@sprintf("Increasing timestep to %.2f days", Δt/86400)) =#
            #=     evolutionLHS = lu(getEvolutionLHS(Δt, diffMat, bdyFluxMat, bottomBdy, topBdy)) =#
            #= end =#
        end
        if i % nStepsPlot == 0
            # plot flow
            iImg += 1
            plotCurrentState(t, chi, uξ, uη, uσ, b, iImg)

            #= # save data =#
            #= file = h5open("b.h5", "w") =#
            #= write(file, "b", b) =#
            #= write(file, "t", t) =#
            #= close(file) =#
        end
        if i % nStepsProfile == 0
            profilePlot(ax, uξ, uη, uσ, b, iξ, tDays)
        end
    end

    ax[1, 1].legend()
    savefig("profiles.png")

    b = reshape(bVec, nξ, nσ)

    return b
end
b = evolveFullNL(500)
