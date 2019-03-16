#!/usr/local/bin/julia

using DelimitedFiles

include("input.jl")
include("initialize.jl")
include("integrator.jl")
include("histogram.jl")

# ==============================================================================
#                           Execute simulation
# ==============================================================================
function main()
    X = initialize_box()
    Z = clist()
    
    compute_force!(X, Z)
    kinetic!(X)
    
    Nsteps = trunc(Int, N/R)
    h = R*dt
    
    for i = 1:Nsteps
        EM!(h, X, i, Z)
        kinetic!(X)
        
        if mod(i, outFreq) == 0
            println("Step $i | Max. F = $(maximum(X.f)) | Max. p = $(maximum(X.p)) | PE = $(X.pe) | KE = $(X.ke) | T = $(X.temp)")
        end
    end
    
    histogram(X.sp)
end

# ==============================================================================
@time main()
