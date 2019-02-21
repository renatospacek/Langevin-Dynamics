#!/usr/local/bin/julia
using InteractiveUtils
include("input.jl")
include("initialize.jl")
include("integrator.jl")
include("histogram.jl")

# ==============================================================================
#                           Execute simulation
# ==============================================================================

function main()
#    println("________________________________")
#    println("dt = $(dt) | $(nPart) particles | Number of cells = $(Ncel)")
#    println("________________________________")
    println("Initializing positions and velocities...")
    X = initialize_box()
    
    println("Computing initial forces...")
    X = compute_force(X)
    
    X = kinetic(X)
    
    dW = sqrt(dt)*randn(dim, nPart*N)
    
    println("Starting integration...")
    Nsteps = trunc(Int, N/R)

    for i = 1:Nsteps
        X = EM(R, dW, X, i)
        X = kinetic(X)
        
        if mod(i, outFreq) == 0
            println("Step $i | Max. F = $(maximum(X.f)) | Max. p = $(maximum(X.p)) | PE = $(X.pe) | KE = $(X.ke) | T = $(X.temp) | Tot. E $(X.ke + X.pe)")
        end
    end
    
    println("Plotting histogram...")
    histogram(X.sp)
    
end


function vscale(p::Array{Float64})
    tfac = 3*nPart/β
    fac = sqrt(tfac/sum(sum(p.*p)))
    p *= fac
    
    return p
end


@time main()









