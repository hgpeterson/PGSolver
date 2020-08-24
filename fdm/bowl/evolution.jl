using SparseArrays, PyPlot, LinearAlgebra, Printf

close("all")
plt.style.use("~/presentation_plots.mplstyle")

include("setParams.jl")
include("myJuliaLib.jl")
include("plottingLib.jl")
include("inversion.jl")

"""
    matrices = getEvolutionMatrices()

Compute the matrices needed for evolution equation integration.
"""
function getEvolutionMatrices()
    nPts = nρ*nσ

    umap = reshape(1:nPts, nρ, nσ)    
    diffMat = Tuple{Int64,Int64,Float64}[]         # diffusion operator matrix 
    diffVec = zeros(nPts)                          # diffusion operator vector 
    bdyFluxMat = Tuple{Int64,Int64,Float64}[]      # flux at boundary matrix
    ρDerivativeMat = Tuple{Int64,Int64,Float64}[]  # advection operator matrix (ρ)
    σDerivativeMat = Tuple{Int64,Int64,Float64}[]  # advection operator matrix (σ)

    # Main loop, insert stencil in matrices for each node point
    for i=1:nρ
        # periodic in ρ
        iL = mod1(i-1, nρ)
        iR = mod1(i+1, nρ)

        # interior nodes only for operators
        for j=2:nσ-1
            row = umap[i, j] 

            # dσ stencil
            fd_σ = mkfdstencil(σ[j-1:j+1], σ[j], 1)
            κ_σ = sum(fd_σ.*κ[i, j-1:j+1])

            # dρ stencil
            ii = [i-1 i i+1]
            if i == 1
                ii .+= 1
            elseif i == nρ
                ii .-= 1
            end
            fd_ρ = mkfdstencil(ρ[ii], ρ[i], 1)

            # dσσ stencil
            fd_σσ = mkfdstencil(σ[j-1:j+1], σ[j], 2)

            # diffusion term: dσ(κ(N^2 + dσ(b)/H))/H = N^2*dσ(κ)/H + dσ(κ)*dσ(b)/H^2 + κ*dσσ(b)/H^2
            push!(diffMat, (row, umap[i, j-1], (κ_σ*fd_σ[1] + κ[i, j]*fd_σσ[1])/H(ρ[i])^2))
            push!(diffMat, (row, umap[i, j],   (κ_σ*fd_σ[2] + κ[i, j]*fd_σσ[2])/H(ρ[i])^2))
            push!(diffMat, (row, umap[i, j+1], (κ_σ*fd_σ[3] + κ[i, j]*fd_σσ[3])/H(ρ[i])^2))
            diffVec[row] = N^2*κ_σ/H(ρ[i])

            # ρ advection term: dρ()
            push!(ρDerivativeMat, (row, umap[ii[1], j], fd_ρ[1]))
            push!(ρDerivativeMat, (row, umap[ii[2], j], fd_ρ[2]))
            push!(ρDerivativeMat, (row, umap[ii[3], j], fd_ρ[3]))

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
        push!(bdyFluxMat, (row, umap[i, 1], fd_σ[1]/H(ρ[i])))
        push!(bdyFluxMat, (row, umap[i, 2], fd_σ[2]/H(ρ[i])))
        push!(bdyFluxMat, (row, umap[i, 3], fd_σ[3]/H(ρ[i])))

        # flux at boundaries: top
        row = umap[i, nσ] 
        # dσ stencil
        fd_σ = mkfdstencil(σ[nσ-2:nσ], σ[nσ], 1)
        # flux term: dσ(b)/H = -N^2
        push!(bdyFluxMat, (row, umap[i, nσ-2], fd_σ[1]/H(ρ[i])))
        push!(bdyFluxMat, (row, umap[i, nσ-1], fd_σ[2]/H(ρ[i])))
        push!(bdyFluxMat, (row, umap[i, nσ],   fd_σ[3]/H(ρ[i])))
    end

    # Create CSC sparse matrix from matrix elements
    diffMat = sparse((x->x[1]).(diffMat), (x->x[2]).(diffMat), (x->x[3]).(diffMat), nPts, nPts)
    bdyFluxMat = sparse((x->x[1]).(bdyFluxMat), (x->x[2]).(bdyFluxMat), (x->x[3]).(bdyFluxMat), nPts, nPts)
    ρDerivativeMat = sparse((x->x[1]).(ρDerivativeMat), (x->x[2]).(ρDerivativeMat), (x->x[3]).(ρDerivativeMat), nPts, nPts)
    σDerivativeMat = sparse((x->x[1]).(σDerivativeMat), (x->x[2]).(σDerivativeMat), (x->x[3]).(σDerivativeMat), nPts, nPts)

    return diffMat, diffVec, bdyFluxMat, ρDerivativeMat, σDerivativeMat
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
    nPts = nρ*nσ

    # initial condition
    b = zeros(nρ, nσ)
    chi = zeros(nρ, nσ)
    chiEkman = zeros(nρ, nσ)
    u = zeros(nρ, nσ)
    v = zeros(nρ, nσ)
    w = zeros(nρ, nσ)
    uρ = zeros(nρ, nσ)
    uσ = zeros(nρ, nσ)
    
    # plot initial state of all zeros and no flow
    iImg = 0
    plotCurrentState(0, chi, chiEkman, u, v, w, b, iImg)

    # initialize plot of profiles
    iρ = convert(Int64, round(nρ/2))
    ax = profilePlotInit()

    # flatten for matrix mult
    bVec = reshape(b, nPts, 1)
    #= uVec = reshape(u, nPts, 1) =#
    #= wVec = reshape(w, nPts, 1) =#
    uρVec = reshape(uρ, nPts, 1)
    uσVec = reshape(uσ, nPts, 1)
    umap = reshape(1:nPts, nρ, nσ)    
    bottomBdy = umap[:, 1]
    topBdy = umap[:, nσ]

    # timestep
    Δt = 1
    #= Δt = 10 =#
    nStepsInvert = 10
    nStepsPlot = 10
    nStepsProfile = 100
    adaptiveTimestep = false

    # get matrices and vectors
    diffMat, diffVec, bdyFluxMat, ρDerivativeMat, σDerivativeMat = getEvolutionMatrices()

    # left-hand side for evolution equation (save LU decomposition for speed)
    evolutionLHS = lu(getEvolutionLHS(Δt, diffMat, bdyFluxMat, bottomBdy, topBdy))

    # left-hand side for inversion equations
    inversionLHS = lu(getInversionLHS())

    # vector of H values
    HVec = H.(reshape(r, nPts, 1))

    # main loop
    t = 0
    for i=1:nSteps
        t += Δt

        # implicit euler diffusion
        diffRHS = bVec + diffVec*Δt


        #= # function to compute advection RHS =#
        #= # (note the parentheses here to allow for sparse matrices to work first) =#
        #= fAdvRHS(bVec) = -(uρVec.*(ρDerivativeMat*bVec) + uσVec.*(σDerivativeMat*bVec) + uσVec.*HVec*N^2) =#

        #= # explicit timestep for advection =#
        #= advRHS = explicitRHS(Δt, bVec, fAdvRHS) =#

        # sum the two
        #= evolutionRHS = diffRHS + advRHS =#
        evolutionRHS = diffRHS 

        # no flux boundaries
        evolutionRHS[bottomBdy] .= -N^2
        evolutionRHS[topBdy] .= -N^2

        # solve
        bVec = evolutionLHS\evolutionRHS

        # log
        println(@sprintf("t = %d s (i = %d)", t, i))
        if i % nStepsInvert == 0
            # reshape
            b = reshape(bVec, nρ, nσ)

            # invert buoyancy for flow
            chi, u, v, w, U = invert(b, inversionLHS)
            uρ = u
            uσ = (w - repeat(σ', nρ, 1).*Hr.(r).*u)./H.(r)
            #= uVec = reshape(u, nPts, 1) =#
            #= wVec = reshape(w, nPts, 1) =#
            uρVec = reshape(uρ, nPts, 1)
            uσVec = reshape(uσ, nPts, 1)
            uCFL = minimum(abs.(dx./u))
            wCFL = minimum(abs.(dz./w))
            #= println(@sprintf("CFL u: %.2f days", uCFL/86400)) =#
            #= println(@sprintf("CFL w: %.2f days", wCFL/86400)) =#
            #= if 0.5*minimum([uCFL, wCFL]) < Δt && adaptiveTimestep =#
            #=     # need to have smaller step size by CFL =#
            #=     Δt = 0.5*minimum([uCFL, wCFL]) =#
            #=     println(@sprintf("Decreasing timestep to %.2f days", Δt/86400)) =#
            #=     evolutionLHS = lu(getEvolutionLHS(Δt, diffMat, bdyFluxMat, bottomBdy, topBdy)) =#
            #= elseif 0.5*minimum([uCFL, wCFL]) > 10*Δt && adaptiveTimestep =#
            #=     # could have much larger step size by CFL =#
            #=     Δt = minimum([0.5*minimum([uCFL, wCFL]), 10*86400]) =#
            #=     println(@sprintf("Increasing timestep to %.2f days", Δt/86400)) =#
            #=     evolutionLHS = lu(getEvolutionLHS(Δt, diffMat, bdyFluxMat, bottomBdy, topBdy)) =#
            #= end =#
        end
        if i % nStepsPlot == 0
            # plot flow
            iImg += 1
            br = rDerivativeTF(b)
            chiEkman = getChiEkman(br)
            errorEkman = norm(chi - chiEkman, 2)
            println(@sprintf("chiEkman L2 error = %.2e", errorEkman))
            plotCurrentState(t, chi, chiEkman, u, v, w, b, iImg)
        end
        if i % nStepsProfile == 0
            profilePlot(ax, u, v, w, b, iρ, t)
        end
    end

    ax[1, 1].legend()
    savefig("profiles.png")

    b = reshape(bVec, nρ, nσ)

    return b
end
b = evolveFullNL(500)
