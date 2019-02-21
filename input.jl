#!/usr/local/bin/julia

using Printf
using LinearAlgebra
using Statistics
using Polynomials
using LaTeXStrings

# ==============================================================================
#                           Simulation parameters
# ==============================================================================

const t 	        = 1.0        # Total simulation time                          
const N 	        = 10^4        # Number of timesteps
const dt            = t/N         # Timestep                                 
const dim	        = 3           # Dimensionality of the system                                  
const ndim          = 10          # Number of particles per side                                 
const nPart         = ndim^dim    # Number of particles 
const L	            = 18          # Side length of simulation box                                                      
const γ 	        = 1           # Damping                                 
const β 	        = 1           # Inverse temperature                                 
const σ	            = 1           # Lennard-Jones, V(σ) = 0                                
const ε	            = 1           # Lennard-Jones, depth of potential well                                
const rc	        = 2^(1/6)*σ   # Cutoff distance                          
const R             = 1           # Timestep scalar                         
#const p             = 3                                    
#const avgM          = 10     
#const sFreq         = 1
const outFreq       = 100         # Output every outFreq timesteps

# ======= Cell list parameters =======
const M = floor(Int, L/rc)        # Number of cells per side                                              
const rn = L/M                    # Side length of cell                                                                                                                    
const Ncel = M^3                  # Number of cells

# ==============================================================================
#                         Base/Type definitions
# ==============================================================================

Base.show(io::IO, f::Float64) = @printf(io, "%1.3f", f)

mutable struct particle{T<:AbstractFloat}
	q::Array{T}
	p::Array{T}
	f::Array{T}
	sp::Array{T}
	pe::T
	ke::T
	temp::T
end
