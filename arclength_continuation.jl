using LinearAlgebra
using DifferentialEquations

function x(xₚₒ,T,F)
	#integrates system during T
	#F is a vector field for an arbitrary systems of ODEs : ẋ = F(x,t,M)
	time_span = (0.0,T)
	M = []
	problem = ODEProblem(F,xₚₒ,time_span,M)
	solution = solve(problem)
	
	return last(solution)
end 

function H(xₚₒ,T,F)
	#Shooting function
	return x(xₚₒ,T,F) - xₚₒ
end	

function reduce_period(xₚₒ,T,F,tolerance)
	#minds minimum period of a periodic solution
	T = abs(T)
	if T < tolerance #xₚₒ is a fixed point
		return T
	else
		primes = [2,3,5,7,11,13]
		for p = primes
			if norm(H(xₚₒ,T/p,F)) < 5*tolerance
				T /= p
			end 
		end
	end
	return T
end


function jac_shooting(xₒ,T,F :: Function,delta)
	#evaluates jacobian of shooting function at (xₒ,T) using value of vector field F for period derivative and finite-difference aproximation for initial conditions
	#delta is the perturbation to each coordinate
	xT = x(xₒ,T,F)
	M = []
	n = size(xₒ)[1]
	del_T = F(xT,T,M) #vector field at evolved position		
	del_xₒ = zeros(n,n) #derivative matrix of H with respect to xₚₒ	
	for i = 1:n
		delta_vec = vec(zeros(n,1))
		delta_vec[i] = delta #used to perturb ith coordinate
		del_xₒ[:,i] = (x(xₒ+delta_vec,T,F) - xT)/delta
	end
	jac = [del_xₒ-I del_T]; # n x (n+1) jacobian of H
	return jac
end  


function shooting_method(x₀,T,F :: Function, delta, tolerance)
	#finds zeros of the shooting function H() of a system F of ODEs using a multidimensional Newton-Raphson scheme
	#tolerance is the allowed deviation from zero in the shooting residue's norm
	#Newton-Raphson method is supplemented with orthogonality condition
	M = []
	n = size(x₀)[1] #num dimensions
	xₒ = copy(x₀)
	T = copy(T)
	H_T = H(xₒ,T,F)
	i = 0
	while ((norm(H_T) > tolerance) && (i<1000))
		jac_H = jac_shooting(xₒ,T,F,delta)
		xT = x(xₒ,T,F)
		del_T = F(xT,T,M)	 
		A = [jac_H ; [transpose(del_T) 0]] #matrix used to solve for Δ, contains orthogonality condition
		Δ = inv(A)*[-H_T ; 0] #Δ contains corrections to T and zₚₒ
		xₒ += Δ[1:n]
		T += Δ[n+1]
		H_T = H(xₒ,T,F)
		i += 1
	end
	xₚₒ = xₒ
	Tₚ = reduce_period(xₚₒ,T,F,tolerance)
	return (xₚₒ,Tₚ)
end

function shooting_method_modified(x₀,T,F :: Function,delta,tolerance,p)
	#finds zeros of the shooting function H() of a system F of ODEs using a multidimensional Newton-Raphson scheme
	#tolerance is the allowed deviation from zero in the shooting residue's norm
	#Newton-Raphson method is supplemented with orthogonality condition as well as orthogonality to the predictor step p
	M = []
	n = size(x₀)[1] #num dimensions
	xₒ = copy(x₀)
	T = copy(T)
	H_T = H(xₒ,T,F)
	i = 0
	while ((norm(H_T)/norm(xₒ) > tolerance) && (i<1000)) #relative error
		jac_H = jac_shooting(xₒ,T,F,delta)
		xT = x(xₒ,T,F)
		del_T = F(xT,T,M)	 
		C = [jac_H ; [transpose(del_T) 0]; transpose(p)] #matrix used to solve for Δ, contains orthogolity to both del_T and the predictor step
		Δ = pinv(C)*[-H_T ; 0; 0] #corrections to T and zₚₒ calculated via Moore-Penrose pseudo-inverse 
		xₒ += Δ[1:n]
		T += Δ[n+1]
		H_T = H(xₒ,T,F)
		i += 1
	end
	xₚₒ = xₒ
	Tₚ = reduce_period(xₚₒ,T,F,tolerance)
	return (xₚₒ,Tₚ)
end

function arclength_continuation(xₒ,T,F :: Function,delta,tolerance,s,steps)
	#Continues a family of zeros to the shooting function parametrized by pseudo-arclength 
	n = size(xₒ)[1] #num dimensions
	M = []
	sols = [[copy(xₒ),copy(T)]]
	last_p = [zeros(n);1]
	for j = 1:steps
		jac_H = jac_shooting(xₒ,T,F,delta)
		B  = [jac_H ; ones(1,n+1)]
		p = inv(B)*[zeros(n,1) ; 1.0] #p is orthogonal to every vector in the jacobian
			
		if norm(p) != 0 
			p = s*p/norm(p)
		end
		xₒ += s*p[1:n]
		T  += s*p[n+1]	
         	H_T = H(xₒ,T,F)
		i = 0
		xₚₒ,Tₚ = shooting_method_modified(xₒ,T,F,delta,tolerance,p)
		push!(sols,[copy(xₚₒ),copy(Tₚ)])
		last_p = copy(p)
	end
	return sols
end
			
	

