using LinearAlgebra
using DifferentialEquations

function z(zₚₒ,T,F)
	#integrates system during T
	#F is a vector field for an arbitrary systems of ODEs : ẋ = F(x,t,M)
	time_span = (0.0,T)
	M = []
	problem = ODEProblem(F,zₚₒ,time_span,M)
	solution = solve(problem)
	
	return last(solution)
end 

function H(zₚₒ,T,F)
	#Shooting function
	return z(zₚₒ,T,F) - zₚₒ
end	

function reduce_period(F,zₚₒ,T,tolerance)
	if T < tolerance #zₚₒ is a fixed point
		return T
	else
		primes = [2,3,5,7,11,13]
		for p = primes
			while norm(H(zₚₒ,T/p,F)) < tolerance
				T /= p
			end 
		end
		return T
	end
end

function shooting_method(F :: Function,z₀,T₀,tolerance,delta)
	#finds zeros of the shooting function H() of a system F of ODEs using a multidimensional Newton-Raphson scheme
	#tolerance is the allowed from zero in the norm of shooting residue;
	#delta is the perturbation of each coordinate used for a finite-difference analysis 
	M = []
	n = size(z₀)[1] #num dimensions
	zₚₒ = z₀
	T = T₀
	H_T = H(zₚₒ,T,F)
	i = 0
	while ((norm(H_T) > tolerance) && (i<1000))
		zT = z(zₚₒ,T,F)
		del_T = F(zT,T,M) #vector field at evolved position
		del_zₚₒ = zeros(n,n) #derivative of z with respect to zₚₒ	
		for i = 1:n
			delta_vec = vec(zeros(n,1))
			delta_vec[i] = delta #used to perturb ith coordinate
			del_zₚₒ[:,i] = (z(zₚₒ+delta_vec,T,F) - zT)/delta
		end
		A = [[del_zₚₒ-I del_T]; [transpose(del_T) 0]] #matrix used to solve for Δ
		Δ = -inv(A)*[H_T ; 0] #Δ contains corrections to T and zₚₒ
		zₚₒ += Δ[1:n]
		T += Δ[n+1]
		H_T = H(zₚₒ,T,F)
		i += 1
		print(i," ")
	end	
	T = reduce_period(F,zₚₒ,T,tolerance)
	return (zₚₒ,T)
end
