function F_duffing_hardening2(x, t, M)
	#1DOF system with cubic stiffness with eps = 0.2
	x1, x2 = x #initial conditions
	dx1 = x2
	dx2 = -x1 - 0.2*x1^3
	return [dx1;dx2]
end

function duffing_energy_hardening2(x)
	#conserved energy for 1DOF system with cubic stiffness, gauged so that the minimum is 0
	x1, x2 = x
	return 1/2 * x2*x2 + 1/4 * 0.2*x1*x1*x1*x1 + 1/2* x1*x1
end
