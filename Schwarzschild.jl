#= 
Author: Ksymena Poradzisz
Contact: ksymena.poradzisz@gmail.com
Affiliation: Jagiellonian University
Created: [2024-06-28]
Updated: [2024-08-31]
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
    println("An error occured while importing LIBSchw.jl: $(e)")
end
println("Pass β value:")
β_str = readline()
β = parse(Int, β_str)
println("Pass v value:")
v_str = readline()
v  = parse(Float64, v_str)
println("Give me the dimension of the square you want to obtain in the visualisation: ")
a_str = readline()
a_box = parse(Int, a_str)
#const v = -0.3 # velocity in c units
#const β = 1  # Thermodynamic β=1/(kT), aka coolness
const γ = 1 / sqrt(1 - v^2) #global usage unnecessary
#const infinity = 150
M = 1; m_0 = 1;

r_box = 2*a_box / sqrt(2) 
println("Wymiary: r = $(r_box), a = $(a_box)")


function create_r_tbl(start, ending, step, M)
    seq1 = collect(start:step:ending)
    result = seq1 .* M#filter(x -> !(x in [xi_hor, xi_ph, xi_mb]), seq) # .* operator multiplies sequence by a number elementwise
    return result
end
# Creating tables for values
J_t_ABS_sch_values = Float64[]
J_r_ABS_sch_values = Float64[]
J_φ_ABS_sch_values = Float64[]
J_t_SCATT_sch_values = Float64[]
J_r_SCATT_sch_values = Float64[]
J_φ_SCATT_sch_values = Float64[]
J_X_ABS_sch_values = Float64[]
J_Y_ABS_sch_values = Float64[]
J_X_SCATT_sch_values = Float64[]
J_Y_SCATT_sch_values = Float64[]
J_X_TOTAL_sch_values = Float64[]
J_Y_TOTAL_sch_values = Float64[]
n_values = Float64[] 
timestamps = String[]
r_values = Int64[]
φ_values = Float64[]
x_values = Float64[]
y_values = Float64[]
φ_table = LinRange(-π, π, 20)
r_table = create_r_tbl(5,r_box, 1, M)
for φ in φ_table
    for ksi in r_table
    	println(ksi)
        timestamp = string(Dates.now())
        push!( r_values,ksi)
        push!(φ_values,φ)
        push!(timestamps, timestamp) #date and time for which the data was produced
        #alfa = 0.0001
        #println("ksi ",ksi," φ = ", φ)
        J_t_ABSsch = J_t_ABS_sch(ksi,φ)
        J_r_ABSsch = J_r_ABS_sch(ksi,φ)
        J_φ_ABSsch = J_φ_ABS_sch(ksi,φ)
        J_t_SCATTsch = J_t_SCATT_sch(ksi, φ)
        J_r_SCATTsch = J_r_SCATT_sch(ksi, φ)
        J_φ_SCATTsch = J_φ_SCATT_sch(ksi, φ)
        push!(J_t_ABS_sch_values, J_t_ABSsch)
        push!(J_r_ABS_sch_values, J_r_ABSsch)
        push!(J_φ_ABS_sch_values, J_φ_ABSsch)
        push!(J_t_SCATT_sch_values, J_t_SCATTsch)
        push!(J_r_SCATT_sch_values, J_r_SCATTsch)
        push!(J_φ_SCATT_sch_values, J_φ_SCATTsch)
        x = ksi*cos(φ)
        y = ksi*sin(φ)
        push!(x_values, x)
        push!(y_values, y)
        J_r_TOTALsch = J_r_ABSsch+J_r_SCATTsch
        J_φ_TOTALsch = J_φ_ABSsch+J_φ_SCATTsch
        J_t_TOTALsch = J_t_ABSsch+ J_t_SCATTsch
        J_X_TOTALsch = J_r_TOTALsch * cos(φ) - J_φ_TOTALsch * sin(φ) / (M*ksi) #J^x  = J^r Cosφ - (J^φ) Sinφ /r (bo wskazniki)
        J_Y_TOTALsch = J_r_TOTALsch * sin(φ) + J_φ_TOTALsch  * cos(φ) /(M*ksi) #J^y  = J^r Sinφ + J^φ Cosφ /r (bo wskaznikii)
        push!(J_X_TOTAL_sch_values, J_X_TOTALsch)
        push!(J_Y_TOTAL_sch_values, J_Y_TOTALsch)

        #calculating n - the surface density
        n_s = sqrt(-J_φ_TOTALsch^2/ksi^2 + J_t_TOTALsch^2 *ksi/(-2*M+ksi)-(-2*M + ksi)/ksi * J_r_TOTALsch^2)
        #push!(n_s)
        _α_ = 1
        n_infty = 2*π * _α_ * m_0^3 *(1+β)/β^2 * exp(-β)
        n = n_s/n_infty
        push!(n_values, n)
    end
end
data = DataFrame(r = r_values,φ = φ_values,x = x_values, y = y_values,timestamp = timestamps,
                J_t_ABSsch = J_t_ABS_sch_values, J_r_ABSsch = J_r_ABS_sch_values,J_φ_ABSsch = J_φ_ABS_sch_values,
                J_t_SCATTsch = J_t_SCATT_sch_values, J_r_SCATTsch = J_r_SCATT_sch_values,J_φ_SCATTsch = J_φ_SCATT_sch_values,
                J_X_TOTAlsch = J_X_TOTAL_sch_values, J_Y_TOTALsch = J_Y_TOTAL_sch_values, n = n_values)
#saving data to a file
timestamp_for_file = Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")
filename = "/home/korizekori/magisterka/Schwarzschild/data_Schwarzschild_beta_$(β)_v_$(v)_dim_$(a_box)_$(timestamp_for_file).csv"
if isfile(filename)
    CSV.write(filename, data, append = true)
else
    CSV.write(filename, data)
end
