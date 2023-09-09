using Plots
using LinearAlgebra
using DifferentialEquations

function plot_2d(xₚ::Vector{Float64},T::Float64, F ; dim1::Int64=1,dim2::Int64=2,color="blue",x_lim::Float64=5.0,y_lim::Float64=5.0)
	#integrates system during T and plots result in phase space
	time_span = (0.0,T)
	M = []
	x0, y0 = xₚ[dim1], xₚ[dim2]
	problem = ODEProblem(F,xₚ,time_span,M)
	solution = solve(problem)
	plt = plot(solution,idxs=(dim1,dim2),plotdensity=10000,label="solution",xaxis = "x",yaxis="ẋ",xlim=(-x_lim,x_lim),ylim=(-y_lim,y_lim))
	scatter!([x0],[y0],color=color,label="xₚ") #plots initial condition
	display(plt)
end  

function plot_2d_many(initial_conditions,periods,F ; dim1::Int64=1,dim2::Int64=2,color="black",x_lim::Float64=5.0,y_lim::Float64=5.0)
	#integrates multiple systems during T and plots result in phase space
	n = size(periods)[1]
	T = periods[1]
	time_span = (0.0,T)
	M = []	
	x0, y0 = initial_conditions[1][dim1], initial_conditions[1][dim2]
	problem = ODEProblem(F,initial_conditions[1],time_span,M)
	solution = solve(problem)
	plt = plot(solution,idxs=(dim1,dim2),plotdensity=10000,label="T = $(round(T;digits=2))",xaxis = "x",yaxis="ẋ",xlim=(-x_lim,x_lim),ylim=(-y_lim,y_lim),legend=:outertopright,color=:red)
	scatter!([x0],[y0],color=color,label="") #plots initial condition
	for k = 2:n
		T = periods[k]
		time_span = (0.0,T)
		M = []	
		x0, y0 = initial_conditions[k][dim1], initial_conditions[k][dim2]
		problem = ODEProblem(F,initial_conditions[k],time_span,M)
		solution = solve(problem)
		plot!(plt,solution,idxs=(dim1,dim2),plotdensity=10000,label="T = $(round(T,digits=2))",xaxis = "x",yaxis="ẋ",xlim=(-x_lim,x_lim),ylim=(-y_lim,y_lim),color=:red)
		scatter!([x0],[y0],color=color,label="")
	end
	display(plt)
end  
