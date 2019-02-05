#!/usr/local/bin/julia

using LinearAlgebra
using Statistics
using Polynomials

include("input.jl")
include("algorithm.jl")
include("cell_lists.jl")

# ==============================================================================
#                           Execute simulation
# ==============================================================================
function main()

    X = initialize(dim, nPart, L)                                               # Initialize positions
    
    #head, list, rn = initialize_list(X.q)                                      # Construct neighbor list
                         
    c = strongconv(X)
    
    println("\n \n Strong convergence = $(c[2]) \n \n")
end

@time main()









