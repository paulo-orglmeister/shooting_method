using Plots
using LaTeXStrings
using DifferentialEquations

function plot_4d(x0,T,F;name="my_4d_plot",save_fig=false,ignore_x4=false)
    #plots solution to 4-dimensional ODE system using xyz and color gradient
    timespan = (0.0,T)
    problem = ODEProblem(F,x0,timespan,[])
    solution = solve(problem,dtmax=0.05) #overrides solver dt for better drawings
    x1 = solution[1,:]
    x2 = solution[2,:]
    x3 = solution[3,:]
    x4 = solution[4,:]
    xmin, ymin, zmin = 1.2*minimum(x1), 1.2*minimum(x2), 1.2*minimum(x3)
    xmax, ymax, zmax = 1.2*maximum(x1), 1.2*maximum(x2), 1.2*maximum(x3)
	if ignore_x4 == false 
    	plt = plot(x1,x2,x3,
		    idxs=(1,2,3),
		    plotdensity=10000,
		    xlim=(xmin,xmax),
		    ylim=(ymin,ymax),
		    zlim=(zmin,zmax),
		    colorbar =:top,
		    color=:jet,
		    framestyle=:box,
		    line_z = x4,
		    linewidth = 2.0,
		    legend=false,
			)
	else 
		plt = plot(x1,x2,x3,
		    idxs=(1,2,3),
		    plotdensity=10000,
		    xlim=(xmin,xmax),
		    ylim=(ymin,ymax),
		    zlim=(zmin,zmax),
		    framestyle=:box,
		    linewidth = 2.0,
		    legend=false,
			color=:jet
		    )
	end
	display(plt)
	if save_fig 
		savefig(name)
	end 
end

function plot_4d_many(initial_conditions,periods,F;name="my_4d_plot",save_fig=false,ignore_x4=false)
	#integrates multiple systems during T and plots result in phase space
	n = size(periods)[1]
	timespan = (0.0,periods[1])
	x0 = initial_conditions[1]
    problem = ODEProblem(F,x0,timespan,[])
    solution = solve(problem,dtmax=0.05) #overrides solver dt for better drawings
    x1 = solution[1,:]
    x2 = solution[2,:]
    x3 = solution[3,:]
    x4 = solution[4,:]
    xmin, ymin, zmin = 5.0*minimum(x1), 4.0*minimum(x2), 10.0*minimum(x3)
    xmax, ymax, zmax = 5.0*maximum(x1), 4.0*maximum(x2), 10.0*maximum(x3)
	if ignore_x4 == false 
    	plt = plot(x3, x2, x1,
		    idxs=(1,2,3),
		    plotdensity=10000,
		    xlim=(xmin,xmax),
		    ylim=(ymin,ymax),
		    zlim=(zmin,zmax),
		    colorbar =:top,
		    color=:jet,
		    framestyle=:box,
		    line_z = x4,
		    linewidth = 2.0,
		    legend=false,
			)
	else 
		plt = plot(x1,x2,x3,
		    idxs=(1,2,3),
		    plotdensity=10000,
		    xlim=(xmin,xmax),
		    ylim=(ymin,ymax),
		    zlim=(zmin,zmax),
		    framestyle=:box,
		    linewidth = 2.0,
		    legend=false,
			color=:jet
		    )
	end
	for k = 2:n
		timespan = (0.0,periods[k])
		x0 = initial_conditions[k]
    	problem = ODEProblem(F,x0,timespan,[])
    	solution = solve(problem,dtmax=0.05) #overrides solver dt for better drawings
    	x1 = solution[1,:]
    	x2 = solution[2,:]
    	x3 = solution[3,:]
    	x4 = solution[4,:]
		xmin, ymin, zmin = 22.0*minimum(x1), 10.0*minimum(x2), 22.0*minimum(x3)
		xmax, ymax, zmax = 22.0*maximum(x1), 10.0*maximum(x2), 22.0*maximum(x3)
		if ignore_x4 == false 
    		plt = plot!(x3,x1,x2,
			    idxs=(1,2,3),
			    plotdensity=10000,
			    colorbar =:top,
			    color=:jet,
			    framestyle=:box,
			    line_z = x4,
			    linewidth = 2.0,
			    legend=false,
				)
		else 
			plt = plot!(x1,x2,x3,
			    idxs=(1,2,3),
			    plotdensity=10000,
			    framestyle=:box,
			    linewidth = 2.0,
			    legend=false,
				color=:jet
			    )
		end
	end
	display(plt)
	if save_fig 
		savefig(name)
	end 
end