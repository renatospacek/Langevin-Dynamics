#!/usr/local/bin/julia

include("input.jl")

# ==============================================================================
#			            Compute energy and forces
# ==============================================================================

function compute_force(forces::Array{Float64},
                       coords::Array{Float64},
                       head::Array{Int},
                       list::Array{Int},
                       rn::Float64)
                       
    forces, pe = build_list(forces, coords, head, list, rn)
    
    return forces, pe
end


function kinetic(L::Int64, vels::Array{Float64})
	ke = 0.0
	
	for i = 1:nPart
		ke += 1/2*dot(vels[i], vels[i])
	end
	
	return ke
end

# ==============================================================================
#		            Lennard-Jones Potential and Force
# ==============================================================================

function LJpe(r)
    pe = 4*ε*((σ/r)^12 - (σ/r)^6)
	
	return pe
end

function LJforce(r)
    force = -(4*ε*(12*(σ./r).^13 - 6*(σ./r).^7))
    
	return force
end
