function F_lure(x,t, M)
	α, β = -1.0, 1.0
	x1, x2, x3 = x
	dx1 = x2
	dx2 = x3
	dx3 = -α * x3 - β * x2 - x1 + (x1)^2
	return [dx1;dx2;dx3]
end
