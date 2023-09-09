using LinearAlgebra
using DifferentialEquations

function x(xₒ,T,F)
	#integrates ẋ = F(x,t,M) from xₒ during T. F is a vector field for an arbitrary systems of ODEs
	M = [] #parameters (useless here)
	time_interval = (0.0,T)
	problem = ODEProblem(F,xₒ,time_interval,M)
	solution = solve(problem)
	return last(solution)
end 

function H(xₒ,T,F)
	#Shooting function
	return x(xₒ,T,F) - xₒ
end	

function get_jac_H(xₒ,T,F :: Function,delta)
	#returns jacobian of a shooting function at (xₒ,T). Derivatives with respect to xₒ and T are determined, respectively, by finite-differences and the value of the vector field F	
	#delta is the perturbation to each coordinate. Works better if less than 1.0e-9
	M = []
	n = size(xₒ)[1]
	xT = x(xₒ,T,F)
	del_T = F(xT,T,M) #vector field at evolved position		
	del_xₒ = zeros(n,n) #derivative matrix of H with respect to xₚ	
	for i = 1:n
		delta_vec = vec(zeros(n,1))
		delta_vec[i] = delta #used to perturb i-th coordinate
		del_xₒ[:,i] = (x(xₒ+delta_vec,T,F) - xT)/delta #finite difference
	end
	jac_H = [del_xₒ-I del_T] # n x (n+1) jacobian of H
	return jac_H
end

function reduce_period(F,xₚ,T,tolerance)
	#minds minimum period of a periodic solution
	if T < tolerance #xₚ is a fixed point
		print("Warning: period T < tolerance. x_p assumed to be a fixed point of ODE system")
		return T
	else
		primes = [2,3,5] # usually enough
		for p = primes
			if norm(H(xₚ,(T/p),F))/norm(xₚ) < tolerance/p #divide by p because if T/p is to small x(xₚ,T/p) might be to close to xₚ
				T /= p
			end 
		end
		return T
	end
end


