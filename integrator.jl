#!/usr/local/bin/julia

include("input.jl")
include("cell_lists.jl")
include("PBC.jl")

# ==============================================================================
#                           Euler-Maruyama
# ==============================================================================
function EM!(h::Float64, X::particle{Float64}, i::Int, Z::clist)
    
    @. X.qtmp = X.q + h.*X.p
	X.p +=  h.*X.f .- γ.*X.p.*h .+ Σ*dt2.*randn(dim, nPart)
	X.q = X.qtmp
	
	PBC!(X.q)
	compute_force!(X, Z)
end

