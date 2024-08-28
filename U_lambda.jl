using Plots
using LaTeXStrings

# Function defining U_λ(ξ)
function U_λ(ξ, λsqr)
    temp = (1 - 2/ξ) * (1 + λsqr/ξ^2)
    return temp
end

# Function to draw multiple U_λ(ξ) for different λ^2
function draw_and_save_multiple_U_λ(λ2_values,filename="plot.png")
    plt = plot() 
    plot!([2, 500], [1, 1], linestyle=:dot, label=L"asymptote")  

    for λ in λ2_values
        U_temp(ξ) = U_λ(ξ, λ)
        plot!(plt, U_temp, 2, 500, label=L"\lambda^2\; =\; %$(λ)")  
        if λ >12
            ξ_max = λ/2*(1-sqrt(1-12/(λ)))
            ξ_min = λ/2*(1+sqrt(1-12/(λ)))
            plot!([ξ_max, ξ_max], [0, 10], linestyle=:dot, label=L"ξ_{max}") 
            plot!([ ξ_min,  ξ_min], [0, 10], linestyle=:dot, label=L" ξ_{min}")


        end
    end

    xlabel!(L"\xi") 
    ylabel!(L"U_\lambda(\xi)")  
    title!(L" U_\lambda(\xi)") 

    savefig(plt, filename)

end

draw_and_save_multiple_U_λ([2, 4, 10, 12], "U_lambda_plot.png")