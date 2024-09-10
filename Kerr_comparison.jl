#= 
Author: Ksymena Poradzisz
Contact: ksymena.poradzisz@gmail.com
Affiliation: Jagiellonian University
Updated: [2024-08-31]
=#

using Symbolics
using QuadGK
using DoubleExponentialFormulas
using PolynomialRoots
using Polynomials
using Dates, CSV
using DataFrames
try
    include("LIBKerr.jl")
catch e
    println("An error occurred while importing LIBKerr.jl: $(e)")
end

const v = -0.5 # velocity in c units
const β = 2.0    # Thermodynamic β=1/(kT), aka coolness
const γ = 1 / sqrt(1 - v^2) # global usage unnecessary

const M = 1.0; m_0 = 1.0;

# Define your xticks and yticks
xticks = 10:-0.5:5
yticks = -π:π/6:π

#Creating tables for values
  J_t_ABS_Kerr_values = Float64[]
  J_r_ABS_Kerr_values = Float64[]
  J_φ_ABS_Kerr_values = Float64[]
J_t_SCATT_Kerr_values = Float64[]
J_r_SCATT_Kerr_values = Float64[]
J_φ_SCATT_Kerr_values = Float64[]
r_values = Float64[]
φ_values = Float64[]

for ksi in xticks
    for φ in yticks
        push!(r_values, ksi)
        push!(φ_values, φ)
        println(ksi," ", φ)
        println("J_t_ABS_Kerr")
        J_t_ABS_Kerr_tmp = J_t_ABS_Kerr(ksi, φ, 0.2, m_0)
        println("J_r_ABS_Kerr")
        J_r_ABS_Kerr_tmp = J_r_ABS_Kerr(ksi, φ, 0.2, m_0)
        println("J_φ_ABS_Kerr")
        J_φ_ABS_Kerr_tmp = J_φ_ABS_Kerr(ksi, φ, 0.2, m_0, M)
        println("J_t_SCATT_Kerr")
        J_t_SCATT_Kerr_tmp = J_t_SCATT_Kerr(ksi, φ, 0.2, m_0)
        println("J_r_SCATT_Kerr")
        J_r_SCATT_Kerr_tmp = J_r_SCATT_Kerr(ksi, φ, 0.2, m_0)
        println("J_φ_SCATT_Kerr")
        J_φ_SCATT_Kerr_tmp = J_φ_SCATT_Kerr(ksi, φ, 0.2, m_0, M)

        push!(J_t_ABS_Kerr_values,     J_t_ABS_Kerr_tmp)
        push!(J_r_ABS_Kerr_values,     J_r_ABS_Kerr_tmp)
        push!(J_φ_ABS_Kerr_values,     J_φ_ABS_Kerr_tmp)
        push!(J_t_SCATT_Kerr_values, J_t_SCATT_Kerr_tmp)
        push!(J_r_SCATT_Kerr_values, J_r_SCATT_Kerr_tmp)
        push!(J_φ_SCATT_Kerr_values, J_φ_SCATT_Kerr_tmp)
    end
end

data = DataFrame(r = r_values, φ = φ_values,
                J_t_ABS_Kerr_tmp = J_t_ABS_Kerr_values, J_r_ABS_Kerr_tmp = J_r_ABS_Kerr_values, J_φ_ABS_Kerr_tmp = J_φ_ABS_Kerr_values,
                J_t_SCATT_Kerr_tmp = J_t_SCATT_Kerr_values, J_r_SCATT_Kerr_tmp = J_r_SCATT_Kerr_values, J_φ_SCATT_Kerr_tmp = J_φ_SCATT_Kerr_values)

# Saving data to a file
timestamp_for_file = Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")
current_path = pwd()
file = "DATA_Kerr_comp_$(timestamp_for_file).csv"
filename = joinpath(current_path, file)

if isfile(filename)
    CSV.write(filename, data, append = true)
else
    CSV.write(filename, data)
end
