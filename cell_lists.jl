#!/usr/local/bin/julia

include("input.jl")


# ==============================================================================
#                 Constructs cell lists and performs
# 	 	               interaction computations
# ==============================================================================

# ==============================================================================
#				   TO DO: 
#		- only works in 3D
#		- assumes a cubic simulation box
# ==============================================================================

# L    = side length of simulation box
# rn   = side length of cell
# rc   = interaction cutoff distance
# M    = number of cells per side
# c    = scalar cell index
# Ncel = number of cells (M^3)

# ==============================================================================
#	                         Initialize list
# ==============================================================================

function initialize_list(q::Array{Float64})
	M = floor(Int, L/rc)                                                        # Number of cells per side
	global rn = L/M                                                             # side length of each cell	                                                        
	Ncel = M^3                                                                  # total number of cells

	head::Array{Int64,1} = zeros(Ncel)                                          # holds the index of the first atom in the c-th cell
	list::Array{Int64,1} = zeros(nPart)                                         # holds the atom index to which the i-th atom points

	for i in 1:nPart
		c = cell_index(q[:,nPart])                                              # compute scalar cell index
		list[i] = copy(head[c])                                                 # link to previous occupant (or empty if you're the first)
		head[c] = i
		
	end	
	
	return head, list, rn
end

# ==============================================================================
#	       		               Clear list
# ==============================================================================

function clear_list()
	for i in 1:N
		list[i] = -1;
	end
		
	for i in 1:Ncel
		head[i] = -1;
	end
end

# ==============================================================================
#	                  Build list/compute interactions
# ==============================================================================


function build_list(forces, q, head, list, rn)
    mc = zeros(dim)
    
	for i in 1:N
		for a in 1:3
			mc[a] = floor(Int, q[a,i]/rn)
		end
	
		for mc1 in (mc[1]-1):(mc[1]+1),
			mc2 in (mc[2]-1):(mc[2]+1),
			mc3 in (mc[3]-1):(mc[3]+1)  		                                # Scan the neighbor cells (including itself) of cell c
			
			c1 = cell_index([mc1,mc2,mc3]) 				                                # Calculate the scalar cell index of the neighbor cell
			j = head[c1] 	                                                    # Scan atom j in cell c1
		
			while j != -1	                                                    
				if i < j                                                        # Avoid double counting of pair (i, j)
					r = PBC(q[:,i], q[:,j])
					force = LJforce(r)
					pe[i] += LJpe(r)
					pe[j] -= LJpe(r)
					forces[:,i] += force
					forces[:,j] -= force
				end
			end
			j = list[j]
		end
	end
	
	return forces, pe
end


# ====== Scalar cell index calculator ==========================================

function cell_index(mc::Array{Int})
	c = 1 + mc[3] + mc[2]*M + mc[1]*M^2

	return c
end

function cell_index(q::Array{Float64})
    c = 1 + floor(Int, q[3]/rn + M/2)
          + floor(Int, q[2]/rn + M/2)*M
          + floor(Int, q[1]/rn + M/2)*M^2
      
    return c
end


