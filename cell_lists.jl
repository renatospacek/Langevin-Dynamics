#!/usr/local/bin/julia

include("input.jl")
include("hamiltonian.jl")
include("PBC.jl")

# ==============================================================================
function build_list!(X::particle{Float64}, Z::clist)
   
    for i in 1:nPart
        vec_index!(X, Z, i)  
        c = 1 + Z.mc[1] + Z.mc[2]*M + Z.mc[3]*M^2
        
        @inbounds Z.list[i] = Z.head[c]
        @inbounds Z.head[c] = i
    end
end

# ==============================================================================
function clear_list!(Z::clist)

    fill!(Z.list, 0)
    fill!(Z.head, 0)

end

# ==============================================================================
#                           Compute interactions
# ==============================================================================
function compute_force!(X::particle{Float64}, Z::clist)
    
    fill!(X.f, 0.0)
    X.pe = 0.0
    
    if nopot == 1
        return X
    end
    
    clear_list!(Z)
    build_list!(X, Z)

    for i in 1:nPart       
        vec_index!(X, Z, i)
        
        for ci1 in (Z.mc[1]-1):(Z.mc[1]+1),
            ci2 in (Z.mc[2]-1):(Z.mc[2]+1),
            ci3 in (Z.mc[3]-1):(Z.mc[3]+1)

            adjc = 1 + cind(ci1) + cind(ci2)*M + cind(ci3)*M*M
            j = Z.head[adjc]

            while j != 0
                if i<j
                    update_force!(X, Z, i, j)
                end
                @inbounds j = Z.list[j]
            end
        end
    end
end

# ==============================================================================
function update_force!(X::particle{Float64}, 
                       Z::clist,
                       i::Int64,
                       j::Int64)

    for l in 1:3
        Z.r[l] = X.q[l,i] - X.q[l,j]
    end
    
    dist = PBC(Z)   
    ff = fLJ(dist)
    X.pe += peLJ(dist)

    for l in 1:3    
        X.f[l,i] += ff*Z.r[l]
        X.f[l,j] -= ff*Z.r[l]
    end  
end

# ==============================================================================
function vec_index!(X::particle{Float64}, Z::clist, i::Int)
    for j in 1:dim
        @inbounds Z.dx[j] = remap(X.q[j,i] + L2)
        @inbounds Z.mc[j] = floor(Int, Z.dx[j]/L*M)   
    end
end

# ==============================================================================
function cind(x::Int64)::Int64

    while x >= M; x -= M; end
    while x < 0; x += M; end
    
    return x
end

