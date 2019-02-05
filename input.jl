#!/usr/local/bin/julia

# ==============================================================================
#                           Simulation parameters
# ==============================================================================

# nPart must be must be n^dim, where n ∊ ℕ

t 	    = 1.0                                                               # Total simulation time
N 	    = 2^10                                                              # Number of timesteps
dt          = t/N                                                               # Timestep size
dim	    = 3                                                                 # Dimensionality of the system
nPart       = 8                                                                 # Number of particles
L	    = 1.0                                                               # Side length of simulation box
ρ           = nPart/L^dim                                                       # Density of particles in simulation box
γ 	    = 1.0                                                               # Damping coefficient
β 	    = 1.0                                                               # Inverse temperature
σ	    = 1.0                                                               # LJ parameter | V(σ) = 0
ε	    = 1.0                                                               # LJ parameter | depth of potential well
rc	    = 2^(1/6)*σ                                                         # Interaction cutoff distance      
R           = 1                                                                 # Timestep size multiplier
p           = 3                                                                 # Number of timestep sizes used
M           = 1000                                                              # Number of dW vectors averaged over




# ==============================================================================
#                           Type definitions
# ==============================================================================

mutable struct particle{T<:AbstractFloat}
	q::Array{T}
	p::Array{T}
#	ke::Array{T}
#	pe::Array{T}
end
