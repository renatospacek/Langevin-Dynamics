#!/usr/local/bin/julia

include("input.jl")
include("cell_lists.jl")
include("PBC.jl")

# ==============================================================================
#                           Euler-Maruyama
# ==============================================================================
function EM(R::Int64, dW::Array{Float64,2}, X::particle, i::Int) 
    h = dt*R
	Winc = _addnoise(R, dW, i)
    
	X.q +=  h*X.p
	X.p +=  h*X.f - γ*X.p*h + sqrt(2*γ/β)*Winc
	
	X.q = PBC(X.q)
	X = compute_force(X)
	
	return X
end

# ==============================================================================
#                     Compute additive noise array
# ==============================================================================
function _addnoise(R::Int, dW::Array{Float64}, j::Int)
    Winc = zeros(dim,nPart)
    
    for i in 1:nPart
        Winc[:,i] = sum(dW[:,(nPart*R*(j-1)+R*(i-1)+1):R*(nPart*(j-1)+i)], dims=2)
    end
    
    return Winc
end















