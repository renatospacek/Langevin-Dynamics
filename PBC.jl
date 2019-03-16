#!/usr/local/bin/julia

include("input.jl")

# ==============================================================================
function PBC!(q::Array{Float64,2})
    for j in 1:dim, i in 1:nPart
            q[j,i] = remap(q[j,i]+L2) - L2
    end
end

# ==============================================================================
function PBC(Z::clist)
    
    for j in 1:dim
        Z.r[j] = remap(Z.r[j]+L2) - L2   
    end
    
    dist = sqrt(Z.r[1]^2 + Z.r[2]^2 + Z.r[3]^2)
	lmul!(1/dist, Z.r)
	
    return dist
end

# ==============================================================================
function remap(x::Float64)::Float64

    while x >= L; x -= L; end
    while x < 0; x += L; end
    
    return x
end

















