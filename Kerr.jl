#= 
Author: Ksymena Poradzisz
Contact: ksymena.poradzisz@gmail.com
Affiliation: Jagiellonian University
Created: [2023-10-21]
Updated: [2024-09-08]
Description:
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
const β = 2    # Thermodynamic β=1/(kT), aka coolness
const γ = 1 / sqrt(1 - v^2) # global usage unnecessary

M = 1; m_0 = 1;

# Define your xticks and yticks
xticks = 10:-0.5:5
yticks = -π:π/6:π

# Creating tables for values
J_t_ABS_kerr_values = Float64[]
J_r_ABS_kerr_values = Float64[]
J_φ_ABS_kerr_values = Float64[]
J_t_SCATT_kerr_values = Float64[]
J_r_SCATT_kerr_values = Float64[]
J_φ_SCATT_kerr_values = Float64[]
r_values = Float64[]
φ_values = Float64[]

for ksi in xticks
    for φ in yticks
        push!(r_values, ksi)
        push!(φ_values, φ)
        J_t_ABSkerr = J_t_ABS_kerr(__jt_integrals__,ksi,φ,0,m_0)
        J_r_ABSkerr = J_r_ABS_kerr(__jr_integrals__,ksi,φ,0,m_0)
        J_φ_ABSkerr = J_φ_ABS_kerr(__jφ_integrals__,ksi,φ,0,m_0,M)

        J_t_SCATTkerr = J_t_SCATT_kerr(__jt_integrals__, ksi, φ, 0, m_0)
        J_r_SCATTkerr = J_r_SCATT_kerr(__jr_integrals__, ksi, φ, 0, m_0)
        J_φ_SCATTkerr = J_φ_SCATT_kerr(__jφ_integrals__, ksi, φ,0, m_0, M)
        push!(J_t_ABS_kerr_values, J_t_ABSkerr)
        push!(J_r_ABS_kerr_values, J_r_ABSkerr)
        push!(J_φ_ABS_kerr_values, J_φ_ABSkerr)
        push!(J_t_SCATT_kerr_values, J_t_SCATTkerr)
        push!(J_r_SCATT_kerr_values, J_r_SCATTkerr)
        push!(J_φ_SCATT_kerr_values, J_φ_SCATTkerr)
    end
end

data = DataFrame(r = r_values, φ = φ_values,
                J_t_ABSkerr = J_t_ABS_kerr_values, J_r_ABSkerr = J_r_ABS_kerr_values, J_φ_ABSkerr = J_φ_ABS_kerr_values,
                J_t_SCATTkerr = J_t_SCATT_kerr_values, J_r_SCATTkerr = J_r_SCATT_kerr_values, J_φ_SCATTkerr = J_φ_SCATT_kerr_values)

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
