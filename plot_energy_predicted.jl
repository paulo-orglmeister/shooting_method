using Plots

function plot_energy_predicted(energy_func,sols; title="",logscale=false,epsilon=0.5)
    #plots backbone curve for hardening Duffing oscillator compared with prediction from multiple timing method (Nayfeh) 
	energies = []
	frequencies = [] #ordinary frequency (not angular)
    frequencies_predicted = []
	for k in eachindex(sols)
		push!(energies,energy_func(sols[k][1]))
		push!(frequencies,1/sols[k][2])
	end
    for E in energies
        push!(frequencies_predicted,(1/2pi)*(1+(3/8)*(2E)*epsilon - (15/256)*(2E)^2*epsilon^2)) #expression from Nayfeh (2000)
    end

	if logscale
		scatter(energies,[frequencies,frequencies_predicted],xaxis = "Energy",yaxis="Frequency (Hz)",title=title,titlefontsize=10,markersize=0.5,xscale=:log10,markershape=[:cross],label=["numerical continuation" "prediction"],legend=:bottomright)
	else
		scatter(energies,[frequencies,frequencies_predicted],xaxis = "Energy",yaxis="Frequency (Hz)",title=title,titlefontsize=10,markersize=0.5,markershape=[:cross],label=["numerical continuation" "prediction"],legend=:bottomright)
	end
	#for log scale use magic argument axis =:log
end

