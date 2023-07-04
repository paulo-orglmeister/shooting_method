using Plots

function plot_energy(energy_func,sols; title="",logscale=false)
	energies = []
	frequencies = [] #ordinary frequency (not angular)
	for k in eachindex(sols)
		push!(energies,energy_func(sols[k][1]))
		push!(frequencies,1/sols[k][2])
	end
	if logscale
		scatter(energies,frequencies,xaxis = "Energy",yaxis="Frequency (Hz)",title=title,titlefontsize=10,markersize=0.5,xscale=:log10)
	else
		scatter(energies,frequencies,xaxis = "Energy",yaxis="Frequency (Hz)",title=title,titlefontsize=10,markersize=0.5)
	end
	#for log scale use magic argument axis =:log
end

	

