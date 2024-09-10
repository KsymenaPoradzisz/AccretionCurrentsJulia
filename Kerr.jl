
#= 
Author: Ksymena Poradzisz
Contact: ksymena.poradzisz@gmail.com
Affiliation: Jagiellonian University
Updated: [2024-09-10]
=#
using Pkg; Pkg.activate(".")
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
println("Pass β value:")
β_str = readline()
β = parse(Float64, β_str)
println("Pass v value:")
v_str = readline()
v  = parse(Float64, v_str)
println("Pass α value:")
α_str = readline()
α  = parse(Float64,α_str)
println("Give me the dimension of the square you want to obtain in the visualisation: ")
a_str = readline()
a_box = parse(Float64, a_str)
const γ = 1 / sqrt(1 - v^2) #global usage unnecessary
M = 1; m_0 = 1;
r_box = 2*a_box / sqrt(2) 
function create_r_tbl(start, ending, step, M)
    seq1 = collect(start:step:ending)
    result = seq1 .* M#filter(x -> !(x in [xi_hor, xi_ph, xi_mb]), seq) # .* operator multiplies sequence by a number elementwise
    return result
end
φ_table = LinRange(-π, π, 20)
r_table = create_r_tbl(5,r_box, 1, M)
J_t_ABS_Kerr_values = Float64[]
J_r_ABS_Kerr_values = Float64[]
J_φ_ABS_Kerr_values = Float64[]
J_t_SCATT_Kerr_values = Float64[]
J_r_SCATT_Kerr_values = Float64[]
J_φ_SCATT_Kerr_values = Float64[]
J_X_ABS_Kerr_values = Float64[]
J_Y_ABS_Kerr_values = Float64[]
J_X_SCATT_Kerr_values = Float64[]
J_Y_SCATT_Kerr_values = Float64[]
J_X_TOTAL_Kerr_values = Float64[]
J_Y_TOTAL_Kerr_values = Float64[]
n_values = Float64[] 
timestamps = String[]
r_values = Int64[]
φ_values = Float64[]
x_values = Float64[]
y_values = Float64[]
for φ in φ_table
    for ksi in r_table
        J_t_ABS_Kerr_value = J_t_ABS_Kerr(ksi, φ, α, m_0)
        J_r_ABS_Kerr_value = J_r_ABS_Kerr(ksi, φ, α, m_0)
        J_φ_ABS_Kerr_value = J_φ_ABS_Kerr(ksi, φ, α, m_0, M)
        J_t_SCATT_Kerr_value = J_t_SCATT_Kerr(ksi, φ, α, m_0)
        J_r_SCATT_Kerr_value = J_r_SCATT_Kerr(ksi, φ, α, m_0)
        J_φ_SCATT_Kerr_value = J_φ_SCATT_Kerr(ksi, φ, α, m_0, M)
        
        # Push the evaluated values into the respective arrays
        push!(J_t_ABS_Kerr_values, J_t_ABS_Kerr_value)
        push!(J_r_ABS_Kerr_values, J_r_ABS_Kerr_value)
        push!(J_φ_ABS_Kerr_values, J_φ_ABS_Kerr_value)
        push!(J_t_SCATT_Kerr_values, J_t_SCATT_Kerr_value)
        push!(J_r_SCATT_Kerr_values, J_r_SCATT_Kerr_value)
        push!(J_φ_SCATT_Kerr_values, J_φ_SCATT_Kerr_value)
        
        # Calculate x and y coordinates
        x = ksi * cos(φ)
        y = ksi * sin(φ)
        push!(x_values, x)
        push!(y_values, y)
        
        # Calculate total currents
        J_r_TOTAL_Kerr = J_r_ABS_Kerr_value + J_r_SCATT_Kerr_value
        J_φ_TOTAL_Kerr = J_φ_ABS_Kerr_value + J_φ_SCATT_Kerr_value
        J_t_TOTAL_Kerr = J_t_ABS_Kerr_value + J_t_SCATT_Kerr_value
        
        # Calculate the components of J^x and J^y
        J_X_TOTAL_Kerr = J_r_TOTAL_Kerr * cos(φ) - J_φ_TOTAL_Kerr * sin(φ) / (M * ksi)  # J^x  = J^r Cosφ - (J^φ) Sinφ /r
        J_Y_TOTAL_Kerr = J_r_TOTAL_Kerr * sin(φ) + J_φ_TOTAL_Kerr * cos(φ) / (M * ksi)  # J^y  = J^r Sinφ + J^φ Cosφ /r
        push!(J_X_TOTAL_Kerr_values, J_X_TOTAL_Kerr)
        push!(J_Y_TOTAL_Kerr_values, J_Y_TOTAL_Kerr)
        
        # Calculate n - the surface density
        n_s = sqrt(-J_φ_TOTAL_Kerr^2 / ksi^2 + J_t_TOTAL_Kerr^2 * ksi / (-2 * M + ksi) - (-2 * M + ksi) / ksi * J_r_TOTAL_Kerr^2)
        _α_ = 1
        n_infty = 2 * π * _α_ * m_0^3 * (1 + β) / β^2 * exp(-β)
        n = n_s / n_infty
        push!(n_values, n)
    end
end
    
    # Create DataFrame with the computed values
    data = DataFrame(r = r_values, φ = φ_values, x = x_values, y = y_values, timestamp = timestamps,
                     J_t_ABS_Kerr = J_t_ABS_Kerr_values, J_r_ABS_Kerr = J_r_ABS_Kerr_values, J_φ_ABS_Kerr = J_φ_ABS_Kerr_values,
                     J_t_SCATT_Kerr = J_t_SCATT_Kerr_values, J_r_SCATT_Kerr = J_r_SCATT_Kerr_values, J_φ_SCATT_Kerr = J_φ_SCATT_Kerr_values,
                     J_X_TOTAL_Kerr = J_X_TOTAL_Kerr_values, J_Y_TOTAL_Kerr = J_Y_TOTAL_Kerr_values, n = n_values)
# Saving data to a file
timestamp_for_file = Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")
current_path = pwd()
file = "data_Kerr_beta_$(β)_v_$(v)_alfa_$(α)_dim_$(a_box)_$(timestamp_for_file).csv"
filename = joinpath(current_path, file)

