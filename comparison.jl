#=
Author: Ksymena Poradzisz
Updated: [2024-08-30]
=#


using Symbolics
using QuadGK
using DoubleExponentialFormulas
using PolynomialRoots
using Polynomials
using Dates, CSV
using DataFrames
try
    include("LIBSchw.jl")
catch e
    println("An error occured while importing LIBSchw.jl: $(e)")
end
try
    include("LIBKerr.jl")
catch e
    println("An error occured while importing LIBSchw.jl: $(e)")
end

const v = -0.5 # predkosc w jedn. c
const β = 2  # Thermodynamic β=1/(kT), aka coolness
const γ = 1 / sqrt(1 - v^2) #global usage unnecessary

M = 1; m_0 = 1;


# Define your xticks and yticks
xticks = 10:-0.5:5
yticks = -π:π/6:π

# Initialize arrays to store the values
J_t_SCATT_Rel_values = Float64[]
J_r_SCATT_Rel_values = Float64[]
J_φ_SCATT_Rel_values = Float64[]
J_t_ABS_Rel_values = Float64[]
J_r_ABS_Rel_values = Float64[]
J_φ_ABS_Rel_values = Float64[]
φ_values = Float64[]
r_values = Float64[]  

# Compute the values
for x in xticks
    for ϕ in yticks
        try
            push!(r_values, x)
            push!(φ_values, ϕ)

            J_t_ABS_Rel = J_t_ABS_kerr(__jt_integrals__,x,ϕ,0,1)   / J_t_ABS_sch(x,ϕ)
            J_r_ABS_Rel = J_r_ABS_kerr(__jr_integrals__,x,ϕ,0,1)   / J_r_ABS_sch(x,ϕ)
            J_φ_ABS_Rel = J_φ_ABS_kerr(__jφ_integrals__,x,ϕ,0,1,1) / J_φ_ABS_sch(x,ϕ)

            J_t_SCATT_Rel = J_t_SCATT_kerr(__jt_integrals__, x, ϕ, 0, 1) / J_t_SCATT_sch(x, ϕ)
            J_r_SCATT_Rel = J_r_SCATT_kerr(__jr_integrals__, x, ϕ, 0, 1) / J_r_SCATT_sch(x, ϕ)
            J_φ_SCATT_Rel = J_φ_SCATT_kerr(__jφ_integrals__, x, ϕ,0, 1, 1)/ J_φ_SCATT_sch(x, ϕ)
            push!(J_t_ABS_Rel_values,J_t_ABS_Rel )
            push!(J_r_ABS_Rel_values,J_r_ABS_Rel )
            push!(J_φ_ABS_Rel_values,J_φ_ABS_Rel )
            push!(J_t_SCATT_Rel_values,J_t_SCATT_Rel )
            push!(J_r_SCATT_Rel_values,J_r_SCATT_Rel )
            push!(J_φ_SCATT_Rel_values,J_φ_SCATT_Rel )

        catch e
            println("An error occurred: $e for r = $x and phi = $ϕ")
            push!(J_t_ABS_Rel_values,0 )
            push!(J_r_ABS_Rel_values,0 )
            push!(J_φ_ABS_Rel_values,0 )
            push!(J_t_SCATT_Rel_values,0)
            push!(J_r_SCATT_Rel_values,0)
            push!(J_φ_SCATT_Rel_values,0)
        end
    end
end
data = DataFrame(r = r_values, φ = φ_values,J_t_ABS_Rel = J_t_ABS_Rel_values, J_r_ABS_Rel = J_r_ABS_Rel_values, J_φ_ABS_Rel = J_φ_ABS_Rel_values,J_t_SCATT_Rel = J_t_SCATT_Rel_values, J_r_SCATT_Rel = J_r_SCATT_Rel_values, J_φ_SCATT_Rel = J_φ_SCATT_Rel_values)
#saving data to a file
timestamp_for_file = Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")
filename = "comparison_data_$(timestamp_for_file).csv"
if isfile(filename)
    CSV.write(filename, data, append = true)
else
    CSV.write(filename, data)
end
