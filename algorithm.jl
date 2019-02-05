#!/usr/local/bin/julia

include("plot.jl")
include("input.jl")
include("hamiltonian.jl")

# ==============================================================================
# 				Defines the following functions:
#		initialize() - initializes box and particles
#       weakconv()   - computes weak convergence of integrator
#       strongconv() - computes strong convergence of integrator
#       PBC()        - applies periodic boundary conditions
#		run()		 - runs simulation
# ==============================================================================

# ==============================================================================
#                        TO DO:
#                          - LJforce in EM()           
#                          - Fix initialize() for dim<3        
#                          - cell_list issue in LJforce
#                          - compute energy
#                          - scale velocity with temperature
#                          - additive noise optimization
#                          - PBC?
# ==============================================================================


# ==============================================================================
#                       Initialize positions
# ==============================================================================

function initialize(dim::Int, nPart::Int, L::Float64)
	ndim = trunc(Int, nPart^(1/dim))                                            # Number of particles per side
	dx = L/(ndim)                                                               # Initial distance between particles
	
	q::Array{Float64} = zeros(dim, nPart)
    p = _initialize_vel()
    
	X = particle(q, p)

    l = 1
	for qx in 1:ndim, qy in 1:ndim, qz in 1:ndim                                # Initialize positions
		X.q[:,l] = dx*[qx-1, qy-1, qz-1]
		l += 1
	end
	return X
end

# ==============================================================================
#                        Initialize velocities
# ==============================================================================

function _initialize_vel()
    p = randn(dim, nPart)                                                       # Initialize random velocities
	p = 2*p - ones(size(p))
	sumv = zeros(dim)
	
	for i in 1:3                                                                # zero net momentum
        sumv[i] = sum(p[i,1:nPart])
        p[i,1:nPart] -= ones(nPart)*sumv[i]/nPart
    end
    
    return p
end 
    
# ==============================================================================
#                           Weak Convergence
# ==============================================================================
 
function weakconv()
    avg = zeros(p)
    
    for j = 1:M
		dW = sqrt(dt)*randn(dim, nPart)
		
		for l in 1:p
			dt_vals[l] = dt*R*2^(l-1)		
			q_final[l] = EM(R*2^(l-1), dW, X)
			avg[l] += (q_final[l] - avg[l])/j	
		end

        for ll in 2:p
			avgerr[ll-1] = abs(avg[ll-1] - avg[ll])
		end
	end

	x = dt_vals[3:end]	
	y = copy(avgerr[2:end])
	a = range(x[1], length=100, stop=x[end])

	p1 = polyfit(log.(x), log.(y), 1)
	c = coeffs(p1)
end

# ==============================================================================
#                           Strong Convergence
# ==============================================================================

function strongconv(X::particle)
    errm = zeros(p)
    dt_vals = zeros(p)
    q_final = zeros(p)
    err = zeros(p-1)
    errm = zeros(p-1)
    
    for j = 1:M
		dW = sqrt(dt)*randn(dim, nPart*N)
		
		for l in 1:p
			dt_vals[l] = dt*R*2^(l-1)		
			q_final[l] = EM(R*2^(l-1), dW, X)
		end

		for ll in 2:p
			err[ll-1] = abs(q_final[1] - q_final[ll])
			errm[ll-1] += (err[ll-1] - errm[ll-1])/j
		end
	end

	x = dt_vals[2:end]	
	y = copy(errm[1:end])
	a = range(x[1], length=100, stop=x[end])

	p1 = polyfit(log.(x), log.(y), 1)
	c = coeffs(p1)
	
	data_plot(x,y,p1,a,c)
	return c
end

# ==============================================================================
#                           Euler-Maruyama
# ==============================================================================

function EM(R::Int64, dW::Array{Float64,2}, X::particle)
    h = R*dt
    L = trunc(Int, N/R)
    
    forces = zeros(dim, nPart)
    
	for i = 1:L
		Winc = _addnoise(R, dW, L, i)

		#forces, pe = compute_force(forces, X.q)
		
		X.q +=  h*X.p
		X.p +=  h*forces - γ*X.p*h + sqrt(2*γ/β)*Winc
		
		#X.ke[i] = kinetic(L, X.p)
		#X.pe[i] = sum(pe)
	end
	
	return X.q[1,1]
end

# ==============================================================================
#                     Compute additive noise array
# ==============================================================================

function _addnoise(R::Int, dW::Array{Float64}, N::Int, j::Int)
    Winc = zeros(dim,nPart)
    
    for i in 1:nPart
        Winc[:,i] = sum(dW[:,(nPart*R*(j-1)+R*(i-1)+1):R*(nPart*(j-1)+i)], dims=2)
    end
    
    return Winc
end

# ==============================================================================
#                     Periodic Boundary Conditions
# ==============================================================================

function PBC(vec, L)
	vec[1] -= L*round(Int64, vec(1)/L)
	vec[2] -= L*round(Int64, vec(2)/L)
	vec[3] -= L*round(Int64, vec(3)/L)
	
	return vec
end


# ==============================================================================
#                           Run Simulation
# ==============================================================================

function run(R::Int64, M::Int64, p::Int64, X::particle)

    X = EM(R, dW, X)

    return X
end














