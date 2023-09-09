using Plots

function plot_amplitude_predicted(sols; title="",logscale=false,relative_frequency=false,epsilon=0.5,w0=1.0)
    #plots backbone curve (amplitude by ordinary frequency) for hardening Duffing oscillator compared with prediction from multiple timing method (Nayfeh) 
	amplitudes = []
	frequencies = [] #ordinary frequency (not angular)
    frequencies_predicted = []
	for k in eachindex(sols)
		push!(amplitudes,sols[k][1][1]) #assumes phase conditions x[2] = 0
		push!(frequencies,1/sols[k][2])
	end
    for A in amplitudes
        push!(frequencies_predicted,(1/2pi)*(w0+(3/(8*w0))*(A^2)*epsilon - (15/(256*w0^3))*(A^4)*epsilon^2)) #expression from Nayfeh (2000)
    end
	yaxis="Frequência (Hz)"
	if relative_frequency
		frequencies = frequencies/(w0/2pi)
		frequencies_predicted = frequencies_predicted/(w0/2pi)
		yaxis = "Frequência / Frequência natural" 
	end
	if logscale
		scatter(amplitudes,[frequencies,frequencies_predicted],xaxis = "Amplitude máxima (u.a.)",yaxis=yaxis,title=title,titlefontsize=10,markersize=0.5,xscale=:log10,markershape=[:cross],label=["continuação numérica" "múltiplas escalas"],legend=:topleft)
	else
		scatter(amplitudes,[frequencies,frequencies_predicted],xaxis = "Amplitude máxima (u.a.)",yaxis=yaxis,title=title,titlefontsize=10,markersize=0.5,markershape=[:cross],label=["continuação numérica" "múltiplas escalas"],legend=:topleft)
	end
end

