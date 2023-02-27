function F_duffing(x, t, M)
	#1DOF system with cubic stiffness
	x1, x2 = x #initial conditions
	dx1 = x2
	dx2 = 2x1 - x1^3
	return [dx1;dx2]

end 
