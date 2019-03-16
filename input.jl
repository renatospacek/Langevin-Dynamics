#!/usr/local/bin/julia

using Printf
using LinearAlgebra

# ==============================================================================
#                           Simulation parameters
# ==============================================================================
const nopot         = 0
const t 	    = 10.0                          
const N 	    = 10^5    
const dt            = t/N   
const dt2           = sqrt(dt)                             
const dim	    = 3                
const ndim          = 10             
const nPart         = ndim^dim                                 
const γ 	    = 1.0                                         
const β	            = 1.0                  
const Σ             = sqrt(2*γ/β)                         
const σ	            = 1.0           
const ε	            = 1.0                
const rc	    = 2^(1/6)*σ      
const L	            = ndim*rc                              
const R             = 1                  
const outFreq       = 1000    

const tfac = 2/3/nPart
const L2 = L/2

# ======= Cell list parameters =======
const M = floor(Int, L/rc)                      
const rn = L/M                                                                             
const Ncel = M^3

# ==============================================================================
#                         Base/Type definitions
# ==============================================================================
Base.show(io::IO, f::Float64) = @printf(io, "%1.5f", f)

mutable struct particle{T<:AbstractFloat}
	q::Array{T,2}
	qtmp::Array{T,2}
	p::Array{T,2}
	f::Array{T,2}
	sp::Array{T,1}
	pe::T
	ke::T
	temp::T
	function particle(nPart::Int64, β::T) where T<:AbstractFloat
	    q = zeros(T, 3,nPart)
	    qtmp = zeros(T, 3,nPart)
	    p = zeros(T, 3,nPart)
	    f = zeros(T, 3, nPart)
	    sp = Array{T,1}(undef,nPart)
	    pe::T = 0.0
	    ke::T = 0.0
	    temp::T = 1/β  
		
	    new{T}(q,qtmp,p,f,sp,pe,ke,temp)
    end
end

mutable struct clist
    list::Array{Int64,1}
    head::Array{Int64,1}
    mc::Array{Int64,1}
    dr::Array{Float64,1}
    r::Array{Float64,1}
    dx::Array{Float64,1}
    function clist()
        list = Array{Int64,1}(undef, nPart) 
        head = Array{Int64,1}(undef, M^3)
        mc = Array{Int64, 1}(undef,3)
        dr = Array{Float64,1}(undef,3)
        r = Array{Float64,1}(undef,3)
        dx = Array{Float64,1}(undef,3)
        
        new(list,head,mc,dr,r,dx)
    end
end
