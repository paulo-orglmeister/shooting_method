function F_double_well(x, t, M)
	#1DOF system with two potential wells
	x1, x2 = x #initial conditions
	dx1 = x2
	dx2 = 2x1 - x1^3
	return [dx1;dx2]
end 

function double_well_energy(x)
	#conserved energy for 1DOF double-well system 
	x1, x2 = x
	return 1/2 * x2*x2 + 1/4 * x1*x1*x1*x1 - x1*x1
end
