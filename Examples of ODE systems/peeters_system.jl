function F_peeters(x,t,M)
    #2 DOF system with cubic stiffness from Peeters (2009)
    x1, x2, x3, x4 = x #initial conditions
    dx1 = x3
    dx2 = x4
    dx3 = -0.5*x1^3 - 2*x1 + x2
    dx4 = -2*x2 + x1
    return [dx1;dx2;dx3;dx4]
end

function energy_peeters(x)
    x1,x2,x3,x4 = x
    #conserved energy for 2 DOF system with cubic stiffness
    return 0.5*(x3^2 + x4^2) + 0.125*x1^4 + x1^2 + x2^2 - 0.5*(x2^2+x1^2)
end 