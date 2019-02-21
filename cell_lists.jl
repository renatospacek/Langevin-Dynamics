#!/usr/local/bin/julia

include("input.jl")
include("hamiltonian.jl")
include("PBC.jl")

# ==============================================================================
#	                         Initialize list
# ==============================================================================
function initialize_list(X::particle)
    dx = zeros(Float64, 3)
    mc = zeros(Int, 3)
    head = zeros(Int, Ncel)
    list = zeros(Int, nPart)
    
    for i in 1:nPart
        for j in 1:dim
            dx[j] = mod(X.q[j,i] + L/2, L)
            mc[j] = floor(Int, dx[j]/L*M)   
        end
        c =  1 + mc[1] + mc[2]*M + mc[3]*M^2

        list[i] = head[c]
        head[c] = i
    end
    
    return head, list
end

# ==============================================================================
#                           Compute interactions
# ==============================================================================
function compute_force(X::particle)
    pe = 0.0
    X.f = zeros(dim, nPart)
    X.pe = 0.0
    head, list = initialize_list(X)
    mc = zeros(Int, 3)
    dx = zeros(3)
    ci = zeros(Int, 3)

    for i in 1:nPart       
        for j in 1:dim
            dx[j] = mod(X.q[j,i] + L/2, L)
            mc[j] = floor(Int, dx[j]/L*M)
        end
        
        for ci[1] in (mc[1]-1):(mc[1]+1),
            ci[2] in (mc[2]-1):(mc[2]+1),
            ci[3] in (mc[3]-1):(mc[3]+1)

            for j in 1:dim
                ci[j] = mod(ci[j], M)
            end
            
            c = 1 + ci[1] + ci[2]*M + ci[3]*M^2
            j = head[c]

            while j != 0
                if i<j
                    dr = X.q[:,i] - X.q[:,j]
                    r = PBC(dr)   
                    ff, pe = LJ(r) 
                    
                    X.f[:,i] += ff*r
                    X.f[:,j] -= ff*r
                    X.pe += pe
                end
                j = list[j]
            end
        end
    end

    return X
end


