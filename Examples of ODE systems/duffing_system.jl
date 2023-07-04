function F_duffing(x, t, M)
	#1DOF system with cubic stiffness
	x1, x2 = x #initial conditions
	dx1 = x2
	dx2 = 2x1 - x1^3
	return [dx1;dx2]
end 

function duffing_energy(x)
	#conserved energy for 1DOF system with cubic stiffness, gauged so that the minimum is 0
	x1, x2 = x
	return 1+ 1/2 * x2*x2 + 1/4 * x1*x1*x1*x1 - x1*x1
end
