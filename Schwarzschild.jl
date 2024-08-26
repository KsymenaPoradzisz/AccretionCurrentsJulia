#= 
Author: Ksymena Poradzisz
Contact: ksymena.poradzisz@gmail.com
Affiliation: Jagiellonian University
Created: [2024-06-28]
Updated: [2024-08-26]
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

const v = -1/2 # velocity in c units
const β = 2  # Thermodynamic β=1/(kT), aka coolness
const γ = 1 / sqrt(1 - v^2) #global usage unnecessary
const infinity = 150
M = 1; m_0 = 1;

#######SCHWARZSCHILD#############
λ_max_sch(ξ, ε) = ξ * sqrt(ε^2 / (1 - 2 / ξ) - 1)
λ_c_sch(ε) = sqrt(12 / (1 - 4 / (1 + 3 * ε / sqrt(9 * ε^2 - 8))^2))
U_λ_sch(ξ, λ) = (1 - 2 / ξ) * (1 + λ^2 / ξ^2)
function R_sch(ξ, ε, λ)
    temp = ε^2 - U_λ_sch(ξ, λ)
    if temp <= 0
        return Inf
    else
        return temp
    end

end

X_sch(ξ, ε, λ) = quadde(_ξ_-> λ/(_ξ_^2*sqrt(R_sch(_ξ_, ε, λ))),ξ,Inf)[1]

function ε_min_sch(ξ)
    if ξ<= 3
        temp =  Inf
    elseif ξ >= 4
        temp = 1
    elseif ξ <4 && ξ >3
        temp = sqrt((1-2/ξ)*(1+1/(ξ-3)))
    else
        temp=  nothing
    end

    return temp
end

function J_t_ABS_sch(ξ,φ)
    coeff = -2*m_0^3/ξ
  #  coeff = 1
    integral = coeff *  quadgk( ε-> ε* quadgk(λ-> cosh(β*γ*v*sqrt(ε^2-1)*sin(φ)*sin(X_sch(ξ, ε, λ)))*
    exp(-β*γ*(ε+v*sqrt(ε^2-1)*cos(φ)*cos(X_sch(ξ,ε,λ))))/sqrt(R_sch(ξ, ε, λ)),0,λ_c_sch(ε))[1],1,infinity)[1]
    return integral
end
function J_t_SCATT_sch(ξ, φ)
    coeff = -4*m_0^3/ξ
    #coeff = 1
    integral=  coeff*quadgk( 
                ε-> exp(-β*γ*ε) * ε *
                 quadgk( λ ->  1*cosh(β*γ*v*sqrt(ε^2-1)*cos(φ)*cos(X_sch(ξ, ε, λ))) *cosh(β*γ*v*sqrt(ε^2-1)*sin(φ)*sin(X_sch(ξ, ε, λ)))/sqrt(R_sch(ξ, ε, λ)), λ_c_sch(ε), λ_max_sch(ξ,ε))[1],ε_min_sch(ξ), infinity)[1]
    return integral
end
function J_r_ABS_sch(ξ,φ)
    coeff = -2*m_0^3/(ξ-2) 
    #coeff = 1
    integral = coeff *quadgk( ε->  
    quadgk(λ-> cosh(β*γ*v*sqrt(ε^2-1)*sin(φ)*sin(X_sch(ξ, ε, λ)))*exp(-β*γ*(ε+v*sqrt(ε^2-1)*cos(φ)*cos(X_sch(ξ,ε,λ)))),0,λ_c_sch(ε))[1],1,infinity)[1]
    return integral
end
function J_r_SCATT_sch(ξ,φ)
    coeff = 4*m_0^3/(ξ-2)
    integral = coeff * quadgk(ε-> exp(-β*γ*ε) * quadgk(λ-> cosh(β*γ*v*sqrt(ε^2-1)*sin(φ)*sin(X_sch(ξ,ε,λ)))*sinh(β*γ*v*sqrt(ε^2-1)*cos(φ)*cos(X_sch(ξ,ε,λ))),λ_c_sch(ε),λ_max_sch(ξ,ε))[1],ε_min_sch(ξ), infinity)[1]
    return integral
end

function J_φ_ABS_sch(ξ, φ)
    coeff = -2 *M* m_0^3*M / ξ 
   # coeff = 1
    integral = coeff *quadgk(
        ε -> quadgk(
         λ-> λ/sqrt(R_sch(ξ, ε, λ)) * exp(-β*γ*(ε+v*sqrt(ε^2-1)*cos(φ)*cos(X_sch(ξ,ε,λ))))* 
         sinh(β*γ*v*sqrt(ε^2-1)*sin(φ)*sin(X_sch(ξ,ε,λ))),0, λ_c_sch(ε))[1],1,infinity)[1]
    return integral
end
function J_φ_SCATT_sch(ξ, φ)
    coeff = -4 * M * m_0^3*M / ξ
   # coeff = 1
    integral = coeff * quadgk(ε-> exp(-β*γ*ε) *  quadgk(λ-> λ* sinh(β*γ*v*sqrt(ε^2-1)*sin(φ)*sin(X_sch(ξ,ε,λ))) * cosh(β*γ*v*sqrt(ε^2-1)*cos(φ)*cos(X_sch(ξ,ε,λ))) / sqrt(R_sch(ξ, ε, λ)),λ_c_sch(ε),λ_max_sch(ξ,ε))[1], ε_min_sch(ξ), infinity)[1]
    return integral
end


function create_r_tbl(start, ending, step, M)
    seq1 = collect(start:step:ending)
    #seq2 = [2^k for k in 3:10]
    #seq = unique(vcat(seq1, seq2)) #unique - removes duplicates; vcat - adds sequances vertically
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
r_table = create_r_tbl(5,30, 1, M)
for φ in φ_table
    for ksi in r_table
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
filename = "/home/korizekori/magisterka/Schwarzschild/data_Schwarzschild_$(timestamp_for_file).csv"
if isfile(filename)
    CSV.write(filename, data, append = true)
else
    CSV.write(filename, data)
end
