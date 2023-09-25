include("shooting_method.jl")

function get_tangent(xₒ,T,F,delta)
	#Finds tangent vector t to branch of solutions at xₒ, T for system n-dimensional system of ODEs
	n = size(xₒ)[1]
	jac_H = get_jac_H(xₒ,T,F,delta)
	B  = [jac_H ; [zeros(1,n-1) 1 0] ; [1 zeros(1,n-1) 0]] #nth coordinate of x is forced to be 0 (phase condition), while tangent vector's first coordinate is forced to be 1
	t = pinv(B)*[zeros(n,1); 0 ; 1.0] #since t is a null vector of the jacobian, the shooting function doesn't increase in t's direction
	(norm(t) != 0) && (t = t/norm(t)) #tangent vector normalized to have arclength 1
	return t
end

function 	
	#Continues a family of zeros to the shooting function (periodic solutions), taking predictor steps proportional to pseudo-arclength 
	#nth coordinate of xₒ (assumed to be a velocity) must be 0 as a phase condition
	success = true #true if Newton method converged for corrector step converged 
	n = size(xₒ)[1] #num dimensions
	M = []
	sols = [[copy(xₒ),copy(T)]]
	if force_period_increase
		last_t = [zeros(n);1] #period starts increasing
	else 
		last_t = [zeros(n);-1] #period starts decreasing
	end
	j = 1
	verbose && print("Starting pseudo-arclength continuation with s = $s, delta = $delta, $steps steps")
	while (j < steps) && success
		t = get_tangent(xₒ,T,F,delta)		
		(dot(t,last_t) < 0) && (t = -t) #forces tangent vector to keep same direction in solution branch
		xₒ += s*t[1:n]
		T  += s*t[n+1]	
		xₚ,Tₚ,success = shooting_method(xₒ,T,F,delta,tolerance,phase_condition=3,t=t,verbose=verbose)
		if success 
			push!(sols,[copy(xₚ),copy(Tₚ)])
			last_t = copy(t)
			j += 1
		end
	end
	(~success) && print("Error: Modified Newton-Raphson method failed to converge at step j = $j. \n")	
	return sols
end
			
