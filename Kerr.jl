#= 
Author: Ksymena Poradzisz
Contact: ksymena.poradzisz@gmail.com
Affiliation: Jagiellonian University
Created: [2023-10-21]
Updated: [2024-08-30]
Description:
This Julia script is intended to compute integrals for black hole accretion currents J^\mu
=#

#=
To run code:
1. Update: in Linux shell:  juliaup update
2. Install missing packages:
	a) in Linux shell run julia
	b) in Julia REPL (shell) run: 
		import Pkg; Pkg.add("PolynomialRoots")
		import Pkg; Pkg.add("Symbolics")
		import Pkg; Pkg.add("QuadGK")
		import Pkg; Pkg.add("Polynomials")
        import Pkg; Pkg.add("Plots")
3. Run code from Linux shell: julia AccretionIntegrals.jl 

Expected outcome:



=#
using Symbolics
using QuadGK
using DoubleExponentialFormulas
using PolynomialRoots
using Polynomials
using Dates, CSV
using DataFrames
include("LIBKerr.jl")
const v = -1/2 # predkosc w jedn. c
const β = 2  # Thermodynamic β=1/(kT), aka coolness
const γ = 1 / sqrt(1 - v^2) #global usage unnecessary
M = 1; m_0 = 1;


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


for x in 5.5:-0.5:5 # printing differences between two integrals 

    J_t_ABS_kerr__ = J_t_ABS_kerr(__jt_integrals__,x,1,0,1) 
    J_r_ABS_kerr__ = J_r_ABS_kerr(__jr_integrals__,x,1,0,1)  
    J_φ_ABS_kerr__ = J_φ_ABS_kerr(__jφ_integrals__,x,1,0,1,1) 

    println("x = $x , J_t_ABS: Kerr = ", J_t_ABS_kerr__)
    println("x = $x , J_r_ABS: Kerr = ", J_r_ABS_kerr__)
    println("x = $x , J_φ_ABS: Kerr= ", J_φ_ABS_kerr__)    
    println("")

end
