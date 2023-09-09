include("shooting_funcs.jl")

function shooting_method(x₀,T,F :: Function,delta,tolerance;phase_condition=1,t=[],verbose=false)
	#finds zeros of the shooting function H() of a system F of ODEs using a multidimensional Newton-Raphson scheme
	#tolerance is the allowed deviation from zero in the shooting residue's relative norm
	#delta is the perturbation to each coordinate used in a finite difference derivative approximation
	#Newton-Raphson method is supplemented with phase condition (orthogonality or Poincaré) 
	max_steps = 1000 #arbitrary
	M = []
	n = size(x₀)[1] #num dimensions of ODE system
	xₒ ,T = copy(x₀), copy(T)
	H_T = H(xₒ,T,F)
	rel_error = (norm(H_T)/norm(xₒ))
	i = 0	
	while ((rel_error > tolerance) && (i < max_steps))
		jac_H = get_jac_H(xₒ,T,F,delta)	
		if phase_condition == 1 # nth coordinate of xₒ (assumed to be a velocity) is set to be zero
			A = [jac_H; [zeros(1,n-1) 1.0  0.0]] #matrix used to solve for corrections Δ. 
			Δ = -inv(A)*[H_T ; xₒ[n]] #Δ contains corrections to T and xₒ. xₒ[n] is the phase condition, which should be zero 
		elseif phase_condition == 2
			xT = x(xₒ,T,F)
			del_T = F(xT,T,M)
			A = [jac_H; [transpose(del_T) 0.0]] #Poincaré orthogonality condition
			Δ = -inv(A)*[H_T ; del_T]	
		elseif phase_condition == 3 #every step in newton's method the correction is also required to be orthogonal to a tangent vector t 
			C = [jac_H ; [zeros(1,n-1) 1.0  0]; transpose(t)] #matrix used to solve for Δ, contains orthogolity to both del_T and the predictor step (tangent vector)
			Δ = -pinv(C)*[H_T ; xₒ[n]; 0] #uses Moore-Penrose pseudo-inverse because resulting system is overdetermined
		end
		xₒ += Δ[1:n]
		T += Δ[n+1]
		H_T = H(xₒ,T,F)
		rel_error = (norm(H_T)/norm(xₒ)) 
		i += 1
		verbose && print("xₒ = $(round.(xₒ,sigdigits = 5)), T = $(round(T,sigdigits=5)) \n")
	end
	if (rel_error < tolerance)
		verbose && print("Newton-Raphson method converged after i = $i steps within tolerance = $tolerance \n")
		if phase_condition == 3
			xₚ, Tₚ = xₒ, T #bad idea to reduce period during numerical continuation
			return (xₚ,Tₚ,true) #numerical continuation must know whether convergence was achieved
		else
			xₚ, Tₚ = xₒ, reduce_period(F,xₒ,T,tolerance)
			return (xₚ,Tₚ)
		end
	else		
		print("Error: Newton-Raphson method convergence failed for $max_steps steps. \n Last iterate xₒ, T = $(round.(xₒ,sigdigits=5)), $(round(T,sigdigits=5)) \n")
		if phase_condition == 3
			return (xₒ,T,false)
		else
			return xₒ,T
		end
	end
end