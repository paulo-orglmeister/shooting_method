
function F_vanderpol(x,t,M)
	#Van der Pol oscilator
	x1, x2 = x #initial conditions
	dx1 = x2
	dx2 = -x1 + (1-x1^2)*x2
	return [dx1;dx2]
end
