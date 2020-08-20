using SparseArrays, PyPlot, LinearAlgebra, Printf

close("all")
plt.style.use("~/presentation_plots.mplstyle")

include("setParams.jl")
include("myJuliaLib.jl")
include("plottingLib.jl")
include("inversion.jl")

"""
    fx = xDerivativeTF(field)

Compute dx(`field`) in terrian-following coordinates.
Note: dx() = dξ() - dx(H)*σ*dσ()/H
"""
function xDerivativeTF(field)
    # allocate
    fx = zeros(nξ, nσ)

    # dξ(field)
    for j=1:nσ
        # use the fact that ξ is evenly spaced and periodic
        fx[2:end-1, j] = (field[3:end, j] - field[1:end-2, j])/(2*dξ)
        fx[1, j] = (field[2, j] - field[nξ, j])/(2*dξ)
        fx[end, j] = (field[1, j] - field[end-1, j])/(2*dξ)
    end

    # -dx(H)*σ*dσ(field)/H
    for i=1:nξ
        fx[i, :] .-= Hx(ξ[i])*σ.*differentiate(field[i, :], σ)/H(ξ[i])
    end

    return fx
end

"""
    fz = zDerivativeTF(field)

Compute dz(`field`) in terrian-following coordinates.
Note: dz() = dσ()/H
"""
function zDerivativeTF(field)
    # allocate
    fz = zeros(nξ, nσ)

    # dσ(field)/H
    for i=1:nξ
        fz[i, :] .+= differentiate(field[i, :], σ)/H(ξ[i])
    end

    return fz
end

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
    xDerivativeMat = Tuple{Int64,Int64,Float64}[]  # advection operator matrix (x)
    zDerivativeMat = Tuple{Int64,Int64,Float64}[]  # advection operator matrix (z)

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

            # x advection term: dξ(b) - Hx*σ*dσ(b)/H
            push!(xDerivativeMat, (row, umap[i, j-1], -Hx(ξ[i])*σ[j]*fd_σ[1]/H(ξ[i])))
            push!(xDerivativeMat, (row, umap[i, j],   -Hx(ξ[i])*σ[j]*fd_σ[2]/H(ξ[i])))
            push!(xDerivativeMat, (row, umap[i, j+1], -Hx(ξ[i])*σ[j]*fd_σ[3]/H(ξ[i])))
            push!(xDerivativeMat, (row, umap[iR, j],  1.0/(2*dξ)))
            push!(xDerivativeMat, (row, umap[iL, j], -1.0/(2*dξ)))

            # z advection term: dσ(b)/H
            push!(zDerivativeMat, (row, umap[i, j-1], fd_σ[1]/H(ξ[i])))
            push!(zDerivativeMat, (row, umap[i, j],   fd_σ[2]/H(ξ[i])))
            push!(zDerivativeMat, (row, umap[i, j+1], fd_σ[3]/H(ξ[i])))
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
    xDerivativeMat = sparse((x->x[1]).(xDerivativeMat), (x->x[2]).(xDerivativeMat), (x->x[3]).(xDerivativeMat), nPts, nPts)
    zDerivativeMat = sparse((x->x[1]).(zDerivativeMat), (x->x[2]).(zDerivativeMat), (x->x[3]).(zDerivativeMat), nPts, nPts)

    return diffMat, diffVec, bdyFluxMat, xDerivativeMat, zDerivativeMat
end

function eulerRHS(Δt, y, f)
    return Δt*f(y)
end
function midpointRHS(Δt, y, f)
	f1 = f(y)
    f2 = f(y + Δt*f1/2)
	dy = Δt*f2
    return dy
end
function rk4RHS(Δt, y, f)
	f1 = f(y)
    f2 = f(y + Δt*f1/2)
    f3 = f(y + Δt*f2/2)
    f4 = f(y + Δt*f3)
	dy = Δt/6*(f1 + 2*f2 + 2*f3 + f4)
    return dy
end

"""
    plotCurrentState(t, chi, chiEkman, u, v, w, b, iImg)

Plot the buoyancy and velocity state of the model at time `t` using label number `iImg`.
"""
function plotCurrentState(t, chi, chiEkman, u, v, w, b, iImg)
    ridgePlot(chi, b, @sprintf("streamfunction at t = %.1f days", t/86400), L"$\chi$ (m$^2$ s$^{-1}$)")
    savefig(@sprintf("chi%03d.png", iImg))
    close()

    ridgePlot(chiEkman, b, @sprintf("streamfunction theory at t = %.1f days", t/86400), L"$\chi$ (m$^2$ s$^{-1}$)")
    savefig(@sprintf("chiEkman%03d.png", iImg))
    close()

    ridgePlot(b, b, @sprintf("buoyancy perturbation at t = %.1f days", t/86400), L"$b$ (m s$^{-2}$)")
    savefig(@sprintf("b%03d.png", iImg))
    close()

    #= ridgePlot(u, b, @sprintf("cross-ridge velocity at t = %.1f days", t/86400), L"$u$ (m s$^{-1}$)"; vext=5e-5) =#
    ridgePlot(u, b, @sprintf("cross-ridge velocity at t = %.1f days", t/86400), L"$u$ (m s$^{-1}$)")
    savefig(@sprintf("u%03d.png", iImg))
    close()

    #= ridgePlot(v, b, @sprintf("along-ridge velocity at t = %.1f days", t/86400), L"$v$ (m s$^{-1}$)"; vext=2e-2) =#
    ridgePlot(v, b, @sprintf("along-ridge velocity at t = %.1f days", t/86400), L"$v$ (m s$^{-1}$)")
    savefig(@sprintf("v%03d.png", iImg))
    close()

    ridgePlot(w, b, @sprintf("vertical velocity at t = %.1f days", t/86400), L"$w$ (m s$^{-1}$)")
    savefig(@sprintf("w%03d.png", iImg))
    close()
end

#= """ =#
#=     b = evolveDiffusionOnly(nSteps) =#

#= Solve diffusion equation for `b` for `nSteps` time steps. =#
#= """ =#
#= function evolveDiffusionOnly(nSteps) =#
#=     nPts = nξ*nσ =#

#=     # initial condition =#
#=     b = zeros(nξ, nσ) =#
#=     u = zeros(nξ, nσ) =#
#=     v = zeros(nξ, nσ) =#
#=     w = zeros(nξ, nσ) =#
    
#=     # plot initial state of all zeros and no flow =#
#=     iImg = 0 =#
#=     plotCurrentState(0, u, v, w, b, iImg) =#

#=     # flatten for matrix mult =#
#=     bVec = reshape(b, nPts, 1) =#
#=     umap = reshape(1:nPts, nξ, nσ) =#    
#=     bottomBdy = umap[:, 1] =#
#=     topBdy = umap[:, nσ] =#

#=     # timestep =#
#=     Δt = 10*86400 =#

#=     # get matrices and vectors =#
#=     diffMat, diffVec, bdyFluxMat, xDerivativeMat, zDerivativeMat = getEvolutionMatrices() =#

#=     # implicit euler =#
#=     LHS = I/Δt - diffMat =# 

#=     # no flux boundaries =#
#=     LHS[bottomBdy, :] = bdyFluxMat[bottomBdy, :] =#
#=     LHS[topBdy, :] = bdyFluxMat[topBdy, :] =#

#=     # save LU decomposition for speed =#
#=     LHS = lu(LHS) =#

#=     # left-hand side for inversion equations =#
#=     inversionLHS = lu(getInversionLHS()) =#

#=     # main loop =#
#=     for i=1:nSteps =#
#=         # implicit euler diffusion =#
#=         RHS = bVec/Δt + diffVec =#

#=         # no flux boundaries =#
#=         RHS[bottomBdy] .= -N^2 =#
#=         RHS[topBdy] .= -N^2 =#

#=         # solve =#
#=         bVec = LHS\RHS =#

#=         # log =#
#=         println("t = ", i*Δt/86400, " days") =#
#=         if i % 10 == 0 =#
#=             iImg += 1 =#

#=             # reshape =#
#=             b = reshape(bVec, nξ, nσ) =#

#=             # invert buoyancy for flow =#
#=             inversionRHS = getInversionRHS(b) =#
#=             sol = inversionLHS\inversionRHS =#
#=             chi, u, v, w, U = postProcess(sol) =#

#=             # plot flow =#
#=             plotCurrentState(i*Δt, u, v, w, b, iImg) =#
#=         end =#
#=     end =#

#=     b = reshape(bVec, nξ, nσ) =#

#=     return b =#
#= end =#

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
    nPts = nξ*nσ

    # initial condition
    b = zeros(nξ, nσ)
    chi = zeros(nξ, nσ)
    chiEkman = zeros(nξ, nσ)
    u = zeros(nξ, nσ)
    v = zeros(nξ, nσ)
    w = zeros(nξ, nσ)
    
    # plot initial state of all zeros and no flow
    iImg = 0
    plotCurrentState(0, chi, chiEkman, u, v, w, b, iImg)

    # initialize plot of profiles
    iξ = 1
    ax = profilePlotInit()

    # flatten for matrix mult
    bVec = reshape(b, nPts, 1)
    uVec = reshape(u, nPts, 1)
    wVec = reshape(w, nPts, 1)
    umap = reshape(1:nPts, nξ, nσ)    
    bottomBdy = umap[:, 1]
    topBdy = umap[:, nσ]

    # timestep
    #= Δt = 10*86400 =#
    Δt = 10*86400
    #= nSteps100Days = Int64(round(100*86400/Δt)) =#
    #= nSteps1000Days = Int64(round(1000*86400/Δt)) =#

    # get matrices and vectors
    diffMat, diffVec, bdyFluxMat, xDerivativeMat, zDerivativeMat = getEvolutionMatrices()

    # left-hand side for evolution equation (save LU decomposition for speed)
    evolutionLHS = lu(getEvolutionLHS(Δt, diffMat, bdyFluxMat, bottomBdy, topBdy))

    # left-hand side for inversion equations
    inversionLHS = lu(getInversionLHS())

    # main loop
    t = 0
    for i=1:nSteps
        t += Δt
        tDays = t/86400

        # implicit euler diffusion
        diffRHS = bVec + diffVec*Δt


        # function to compute explicit advection 
        # (note the parentheses here to allow for sparse matrices to work first)
        fAdvRHS(bVec) = -(uVec.*(xDerivativeMat*bVec) + wVec.*(zDerivativeMat*bVec) + wVec*N^2)
        #= advRHS = -Δt*(uVec.*(xDerivativeMat*bVec) + wVec.*(zDerivativeMat*bVec) + wVec*N^2) =#

        #= # euler =# 
        #= advRHS = eulerRHS(Δt, bVec, fAdvRHS) =#
        #= # midpoint =# 
        #= advRHS = midpointRHS(Δt, bVec, fAdvRHS) =#
        # RK4
        advRHS = rk4RHS(Δt, bVec, fAdvRHS)

    
        #= # explicit advection (note the parentheses here to allow for sparse matrices to work first) =#
        #= advRHS = -uVec.*(xDerivativeMat*bVec) .- wVec.*(zDerivativeMat*bVec) .- wVec*N^2 =#

        # sum the two
        evolutionRHS = diffRHS + advRHS

        # no flux boundaries
        evolutionRHS[bottomBdy] .= -N^2
        evolutionRHS[topBdy] .= -N^2

        # solve
        bVec = evolutionLHS\evolutionRHS

        # log
        println(@sprintf("t = %.2f days (i = %d)", tDays, i))
        if i % 2 == 0
            # reshape
            b = reshape(bVec, nξ, nσ)

            # invert buoyancy for flow
            chi, u, v, w, U = invert(b, inversionLHS)
            uVec = reshape(u, nPts, 1)
            wVec = reshape(w, nPts, 1)
            uCFL = minimum(abs.(dx./u))
            wCFL = minimum(abs.(dz./w))
            #= println(@sprintf("CFL u: %.2f days", uCFL/86400)) =#
            #= println(@sprintf("CFL w: %.2f days", wCFL/86400)) =#
            if 0.5*minimum([uCFL, wCFL]) < Δt
                # need to have smaller step size by CFL
                Δt = 0.5*minimum([uCFL, wCFL])
                println(@sprintf("Decreasing timestep to %.2f days", Δt/86400))
                evolutionLHS = lu(getEvolutionLHS(Δt, diffMat, bdyFluxMat, bottomBdy, topBdy))
            elseif 0.5*minimum([uCFL, wCFL]) > 10*Δt
                # could have much larger step size by CFL
                Δt = minimum([0.5*minimum([uCFL, wCFL]), 10*86400])
                println(@sprintf("Increasing timestep to %.2f days", Δt/86400))
                evolutionLHS = lu(getEvolutionLHS(Δt, diffMat, bdyFluxMat, bottomBdy, topBdy))
            end
        end
        if i % 10 == 0
            # plot flow
            iImg += 1
            bx = xDerivativeTF(b)
            chiEkman = getChiEkman(bx)
            plotCurrentState(t, chi, chiEkman, u, v, w, b, iImg)
        end
        if i % 100 == 0
            profilePlot(ax, u, v, w, b, iξ, tDays)
        end
    end

    ax[1, 1].legend()
    savefig("profiles.png")

    b = reshape(bVec, nξ, nσ)

    return b
end
b = evolveFullNL(500)
