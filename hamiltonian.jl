#!/usr/local/bin/julia

include("input.jl")

# ==============================================================================
#		   Lennard-Jones Potential and Force
# ==============================================================================
function fLJ(rr::Float64)::Float64
    if rr > rc
        return 0.0
    else
        r1 = σ/rr
        r2 = r1^2
    
        return 48*ε*r2^3*(r2^3 - 1/2)*r1
    end
end

# ==============================================================================
function peLJ(rr::Float64)::Float64
    if rr > rc
        return 0.0
    else
        r1 = σ/rr
        r2 = r1^2
    
        return 4*ε*r2^3*(r2^3 - 1) + 1
    end
end

# ==============================================================================
#                       Kinetic energy computation
# ==============================================================================
function kinetic!(X::particle{Float64})
    p = sum(X.p.*X.p, dims=1)
    p = reshape(p, (length(p)))
    
    X.ke = sum(p/2)
    X.temp = X.ke*tfac
    append!(X.sp, sqrt.(p))
end















