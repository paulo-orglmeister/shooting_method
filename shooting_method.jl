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

function reduce_period(F,xₚₒ,T,tolerance)
	#minds minimum period of a periodic solution
	if T < tolerance #xₚₒ is a fixed point
		return T
	else
		primes = [2,3,5,7,11,13]
		for p = primes
			while norm(H(xₚₒ,T/p,F)) < tolerance
				T /= p
			end 
		end
		return T
	end
end

function shooting_method(x₀,T,F :: Function, tolerance,delta)
	#finds zeros of the shooting function H() of a system F of ODEs using a multidimensional Newton-Raphson scheme
	#tolerance is the allowed deviation from zero in the shooting residue's norm;
	#delta is the perturbation to each coordinate used for a finite-difference approximation
	M = []
	n = size(x₀)[1] #num dimensions
	xₒ = copy(x₀)
	T = copy(T)
	H_T = H(xₒ,T,F)
	i = 0
	while ((norm(H_T) > tolerance) && (i<1000))
		xT = x(xₒ,T,F)
		del_T = F(xT,T,M) #vector field at evolved position
		del_xₒ = zeros(n,n) #derivative matrix of x with respect to xₚₒ	
		for i = 1:n
			delta_vec = vec(zeros(n,1))
			delta_vec[i] = delta #used to perturb ith coordinate
			del_xₒ[:,i] = (x(xₒ+delta_vec,T,F) - xT)/delta
		end
		A = [[del_xₒ-I del_T]; [transpose(del_T) 0]] #matrix used to solve for Δ
		Δ = -inv(A)*[H_T ; 0] #Δ contains corrections to T and zₚₒ
		xₒ += Δ[1:n]
		T += Δ[n+1]
		H_T = H(xₒ,T,F)
		i += 1
	end
	xₚₒ = xₒ
	Tₚ = reduce_period(F,xₚₒ,T,tolerance)
	return (xₚₒ,Tₚ)
end
