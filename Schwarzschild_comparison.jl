#= 
Author: Ksymena Poradzisz
Contact: ksymena.poradzisz@gmail.com
Affiliation: Jagiellonian University
Created: [2024-06-28]
Updated: [2024-09-09]
Description:
This Julia script is intended to compute integrals for Schwarzschild black hole accretion currents J^\mu and save them to *.csv file
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
    error("An error occured while importing LIBSchw.jl: $(e)")
end
const v = -0.5 # predkosc w jedn. c
const β = 2  # Thermodynamic β=1/(kT), aka coolness
const γ = 1 / sqrt(1 - v^2) #global usage unnecessary

M = 1; m_0 = 1;


# Define your xticks and yticks
xticks = 10:-0.5:5
yticks = -π:π/6:π


# Creating tables for values
J_t_ABS_Schw_values = Float64[]
J_r_ABS_Schw_values = Float64[]
J_φ_ABS_Schw_values = Float64[]
J_t_SCATT_Schw_values = Float64[]
J_r_SCATT_Schw_values = Float64[]
J_φ_SCATT_Schw_values = Float64[]
r_values = Float64[]
φ_values = Float64[]

for ksi in xticks
    for φ in yticks
        push!(r_values, ksi)
        push!(φ_values, φ)
        
        J_t_ABS_Schw_value = J_t_ABS_Schw(ksi, φ)
        J_r_ABS_Schw_value = J_r_ABS_Schw(ksi, φ)
        J_φ_ABS_Schw_value = J_φ_ABS_Schw(ksi, φ)
        J_t_SCATT_Schw_value = J_t_SCATT_Schw(ksi, φ)
        J_r_SCATT_Schw_value = J_r_SCATT_Schw(ksi, φ)
        J_φ_SCATT_Schw_value = J_φ_SCATT_Schw(ksi, φ)
        
        push!(J_t_ABS_Schw_values, J_t_ABS_Schw_value)
        push!(J_r_ABS_Schw_values, J_r_ABS_Schw_value)
        push!(J_φ_ABS_Schw_values, J_φ_ABS_Schw_value)
        push!(J_t_SCATT_Schw_values, J_t_SCATT_Schw_value)
        push!(J_r_SCATT_Schw_values, J_r_SCATT_Schw_value)
        push!(J_φ_SCATT_Schw_values, J_φ_SCATT_Schw_value)
    end
end

data = DataFrame(r = r_values, φ = φ_values,
                 J_t_ABS_Schw = J_t_ABS_Schw_values, J_r_ABS_Schw = J_r_ABS_Schw_values, J_φ_ABS_Schw = J_φ_ABS_Schw_values,
                 J_t_SCATT_Schw = J_t_SCATT_Schw_values, J_r_SCATT_Schw = J_r_SCATT_Schw_values, J_φ_SCATT_Schw = J_φ_SCATT_Schw_values)
timestamp_for_file = Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")
current_path = pwd()
file = "DATA_Schw_comp_$(timestamp_for_file).csv"
filename = joinpath(current_path, file)
if isfile(filename)
    CSV.write(filename, data, append = true)
else
    CSV.write(filename, data)
end
