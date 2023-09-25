using Plots
#plotlyjs()

function plot_backbone(sols,energy_func; title="Backbone curve",logscale=false,natural_frequency=1)
	#Graphs backbone curve for periodic solution family given function to calculate energy. Optional logscale for energy 
	energies, frequencies = [], [] #ordinary frequency (not angular) 
	for k in eachindex(sols)
		push!(energies,energy_func(sols[k][1]))
		push!(frequencies,(natural_frequency/sols[k][2]))
	end
	f_axis = "Frequência (Hz)"
	(natural_frequency != 1) && (f_axis = "Frequência / Frequência natural")
	if logscale
		scatter(energies,frequencies,xaxis = "Energia (u.a.)",yaxis=f_axis,title=title,titlefontsize=10,markersize=0.5,xscale=:log10,legend=false)
	else
		scatter(energies,frequencies,xaxis = "Energia (u.a.)",yaxis=f_axis,title=title,titlefontsize=10,markersize=0.5,legend=false)
	end
end

	

