#!/usr/local/bin/julia

include("input.jl")

# ==============================================================================
function PBC(q::Array{Float64,2}) # position between -L/2 and L/2
    for i in 1:nPart
        for j in 1:dim
            q[j,i] = mod(q[j,i]+L/2, L) - L/2   
        end
    end
    
    return q
end

# ==============================================================================
function PBC(q::Array{Float64,1}) # distance between particles
    for j in 1:dim
        q[j] = mod(q[j]+L/2, L) - L/2   
    end
    
    return q
end












