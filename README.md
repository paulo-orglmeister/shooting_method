Implementation of the pseudo-arclength numerical continuation algorithm described in Peeters et. al (2009) [1].
A shooting method to find periodic solutions to non-linear systems is implemented as well. 

To import necessary functions in Julia REPL:

import("parclength_continuation.jl")

Systems of ODEs are implemented in Julia as functions of the coordinates which return the value of a vector field F, as described in DifferentialEquations.jl [2].
The two degree of freedom system used in [1] for instance is:

function F_peeters(x,t,M)
	x1, x2, x3, x4 = x #initial conditions
	dx1 = x3
	dx2 = x4
	dx3 = x2 -0.5*x1^3 - 2*x1
	dx4 = x1 -2*x2
	return [dx1;dx2;dx3;dx4]
end 

To calculate a periodic solutions to this system a initial condition and period, a finite difference and a relative tolerance must be specified:

include("ode_systems/peeters_system.jl")
x0 = [0.04, -0.04, 0.0, 0.0]
T0 = 2pi/sqrt(3)
delta = 1.0e-12
tolerance = 0.0001
xp, Tp = shooting_method(x0,T0,F_peeters,delta,tolerance,verbose=true)

To plot the solution xp, Tp in phase space:

include("ploting_funcs/plot_2d.jl")
plot_2d(xp,Tp,F_peeters)

The pseudo-arclength continuation algorithm is used to find a branch of solutions, with specified number of steps and arlcength increment:

s = 0.01
steps = 200
sols = parclength_continuation(xp,Tp,F_peeters,delta,tolerance,
		force_period_increase = false, verbose=true,s=s,steps=steps) 

It is usually necessary to tune the parameters s and delta to obtain convergence for a greater number of continuation steps. The branch of solutions can be plotted using:

include("ploting_funcs/plot_backbone")
plot_backbone(sols,peeters_energy,logscale=true)



