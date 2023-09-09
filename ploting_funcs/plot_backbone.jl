using Plots

function plot_backbone(sols,energy_func; title="Backbone curve",logscale=false)
	#Graphs backbone curve for periodic solution family given function to calculate energy. Optional logscale for energy 
	energies, frequencies = [], [] #ordinary frequency (not angular) 
	for k in eachindex(sols)
		push!(energies,energy_func(sols[k][1]))
		push!(frequencies,(2pi/sqrt(3))/sols[k][2])
	end
	if logscale
		scatter(energies,frequencies,xaxis = "Energia (u.a.)",yaxis="Frequência / Frequência natural",title=title,titlefontsize=10,markersize=0.5,xscale=:log10,legend=false)
	else
		scatter(energies,frequencies,xaxis = "Energia (u.a.)",yaxis="Frequência (Hz)",title=title,titlefontsize=10,markersize=0.5,legend=false)
	end
end

	

