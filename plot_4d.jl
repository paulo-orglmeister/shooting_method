using Plots
using LaTeXStrings
using DifferentialEquations

function plot_4d(x0,T,F,name="my_4d_plot")
    timespan = (0.0,T)
    problem = ODEProblem(F,x0,timespan,[])
    solution = solve(problem,dtmax=0.05) #overrides solver dt for better drawings
    x1 = solution[1,:]
    x2 = solution[2,:]
    x3 = solution[3,:]
    x4 = solution[4,:]
    xmin, ymin, zmin = 1.2*minimum(x1), 1.2*minimum(x2), 1.2*minimum(x3)
    xmax, ymax, zmax = 1.2*maximum(x1), 1.2*maximum(x2), 1.2*maximum(x3)
    p = plot(x1,x2,x3,
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
	 #xlabel=L"x_1(t)",
	 #ylabel=L"x_2(t)",
	 #zlabel=L"x_3(t)",
         legend=false,
	 #colorbar_title = L"x_4(t)"
	 )
	savefig(name)
    end

#=
opções passadas pelo prof. Franzini
    GuideFont = "Computer Modern" 
    TickFont = "Computer Modern"

    plot(
        x1,
        x2,
        x3,
        #linewidth=1.5mm,
        #right_margin=12,
        #left_margin=2,
        #bottom_margin=3,
        framestyle=:box,
        legend=false,
        xlabel=L"x_1(\tau)",
        ylabel=L"x_2(\tau)",
        zlabel=L"x_3(\tau)",
        colorbar=:best,
        color=:jet,
        line_z=x4, #the x4 variable is used for line color
        #colorbar_title=L"x_4(\tau)",
        #xguidefontsize=GuideFont,
        #xtickfontsize=TickFont,
        #yguidefontsize=GuideFont,
        #ytickfontsize=TickFont,
        #zguidefontsize=GuideFont,
        #ztickfontsize=TickFont,
    )
=#
