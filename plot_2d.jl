using Plots
using LinearAlgebra
using DifferentialEquations

function plot_2d(xₚₒ::Vector{Float64},T::Float64, F ; dim1::Int64=1,dim2::Int64=2,x_zoom::Float64=5.0,y_zoom::Float64=5.0)
	#integrates system during T and plots result in phase space
	time_span = (0.0,T)
	M = []
	x0, y0 = xₚₒ[dim1], xₚₒ[dim2]
	
	problem = ODEProblem(F,xₚₒ,time_span,M)
	solution = solve(problem)
	plot(solution,idxs=(dim1,dim2),plotdensity=10000,label="solution",xaxis = "x",yaxis="ẋ",xlim=(-x_zoom,x_zoom),ylim=(-y_zoom,y_zoom))
	scatter!([x0],[y0],color="blue",label="xₚₒ")
end  
