#!/usr/local/bin/julia

include("input.jl")

# ==============================================================================
#		            Lennard-Jones Potential and Force
# ==============================================================================
function LJ(r::Array{Float64})
    rr = norm(r)
    
    if rr < rc
        r2 = (σ/rr)^2
    
        ff = 48*ε*r2^4*(r2^3 - 1/2)
        pe = 4*ε*r2^3*(r2^3 - 1)
    else
        ff = 0.0
        pe = 0.0
    end
    
	return ff, pe
end

# ==============================================================================
#                       Kinetic energy computation
# ==============================================================================
function kinetic(X::particle)
    p = sum(X.p.*X.p, dims=1)
    p = reshape(p, (length(p)))
    
    X.ke = sum(p/2)
    X.temp = X.ke*2/3/nPart
    append!(X.sp, sqrt.(p))
    
    return X
end










