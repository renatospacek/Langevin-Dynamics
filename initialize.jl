#!/usr/local/bin/julia

include("input.jl")

# ==============================================================================
#                  Generate initial positions and velocities
# ==============================================================================
function initialize_box()                                          
	dx = L/ndim                                                               
	X = particle(nPart, β)
	_initialize_vel!(X)
	kinetic!(X)

	l = 1
	for qx in 0:ndim-1, qy in 0:ndim-1, qz in 0:ndim-1                                
		X.q[1,l] = dx*qx - L2
		X.q[2,l] = dx*qy - L2
		X.q[3,l] = dx*qz - L2
		l += 1
	end

	return X
end

# ==============================================================================
function _initialize_vel!(X::particle)
	X.p = randn(dim, nPart)                                                       
	X.p .= 2*X.p .- ones(size(X.p))
	sumv = Array{Float64,1}(undef,3)
	
	for i in 1:3                                                                
        sumv[i] = sum(X.p[i,1:nPart])
        X.p[i,1:nPart] -= ones(nPart).*sumv[i]./nPart
    end
   
    # Scale velocities according to temperature
    # 1/2 Σ v^2 = 3/2 nPart/β
    tfact = 3*nPart/β
    fac = sqrt(tfact/sum(X.p.*X.p))
    X.p *= fac
end 










