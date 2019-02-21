#!/usr/local/bin/julia

include("input.jl")

# ==============================================================================
#                  Generate initial positions and velocities
# ==============================================================================
function initialize_box()                                          
    dx = L/ndim                                                               
    shift = L/2*ones(3)
                                                             
    q::Array{Float64} = zeros(dim, nPart)
    f::Array{Float64} = zeros(dim, nPart)
    p = _initialize_vel()
    pe = 0.0
    ke = 0.0
    sp = [0.0]
    temp = 1.0
    
    X = particle(q, p, f, sp, pe, ke, temp)
    X = kinetic(X)
    
    l = 1
	for qx in 0:ndim-1, qy in 0:ndim-1, qz in 0:ndim-1                                
		X.q[:,l] = dx*[qx, qy, qz] - shift
		l += 1
	end
	
	return X
end

function _initialize_vel()
    p = randn(dim, nPart)                                                       
	p = 2*p - ones(size(p))
	sumv = zeros(dim)
	
	for i in 1:3                                                                
        sumv[i] = sum(p[i,1:nPart])
        p[i,1:nPart] -= ones(nPart)*sumv[i]/nPart
    end
   
    # Scale velocities according to temperature
    # 1/2 Σ v^2 = 3/2 nPart/β
    tfac = 3*nPart/β
    fac = sqrt(tfac/sum(p.*p))
    p *= fac
    
    return p
end 










