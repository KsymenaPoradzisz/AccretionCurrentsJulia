using Plots
using Colors
using LaTeXStrings

gr()

# Function defining U_λ(ξ), λsqr  = λ^2
function U_λ(ξ, λsqr)
    temp = (1 - 2/ξ) * (1 + λsqr/ξ^2)
    return temp
end

# Function to draw multiple U_λ(ξ) for different λ^2
function draw_U_λ(λ2_values,filename="plot.png")
    plt = plot() 
    n_colors = length(λ2_values)
    seed_colors = [RGB(1, 1, 1), RGB(0, 0, 0)]  

    # generate distinguishable colors
    colors = distinguishable_colors(n_colors, seed_colors, dropseed = true)
   # colors = distinguishable_colors(length(λ2_values),[RGB(1,1,1), RGB(0,0,0)],dropseed = false) 
    plot!([2, 300], [1, 1], linestyle=:dashdotdot, label=L"asymptote", color = "red")  
    tempξ_max,tempξ_min = 0,0
    for (i,λ) in enumerate(λ2_values)
        U_temp(ξ) = U_λ(ξ, λ)
                
        plot!(plt, U_temp, 2, 300, label=L"\lambda^2\; =\; %$(λ)",color=colors[i],legend=:bottomright)  
        line_color = plt.series_list[end][:seriescolor]
        #title!(L" U_\lambda(\xi)") 
        if λ >12
            limit = 10
            ξ_max = λ/2*(1-sqrt(1-12/(λ)))
            ξ_min = λ/2*(1+sqrt(1-12/(λ)))

            plot!([ξ_max, ξ_max], [0, limit], linestyle=:dash, color = line_color,label=L"ξ_{max}") 
            plot!([ ξ_min,  ξ_min], [0, limit], linestyle=:dot,color = line_color, label=L" ξ_{min}")

        end

        ylims!(plt, (0, 1.5))
        lens!([2,25], [0,1.5], inset = (1, bbox(0.22, 0.38, 0.5, 0.5))) #zoom in 
        xlabel!(L"\xi") 
        ylabel!(L"U_\lambda(\xi)")  
    end

  

    savefig(plt, filename)

end

draw_U_λ([13, 16, 25, 20], "U_lambda_plot.png")
