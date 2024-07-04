#= 
Author: Ksymena Poradzisz
Contact: ksymena.poradzisz@gmail.com
Affiliation: Jagiellonian University
Created: [2023-10-21]
Updated: [2024-06-26]
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
using Plots
using PolynomialRoots
using Polynomials
#using DoubleExponentialFormulas, CSV
@variables M, m_0, ξ, ϵ_σ, ϵ_r, ε, α, λ, _ξ_, φ
machine_epsilon = eps() #the smallest number which can be added to 1.0
v = -0.5 # predkosc w jedn. c
β = 2  # Thermodynamic β=1/(kT), aka coolness
#φ = 1   # some "random" polar angle
################FUNCTIONS############

X(ξ, ε, λ, α) = quadgk(_ξ_ -> (λ + (α * ε) / (1 - 2 / _ξ_)) / (((α^2) / (1 - 2 / _ξ_) + _ξ_^2) * sqrt(Complex(ε^2 - (1 - 2 / _ξ_) * (1 + λ^2 / _ξ_^2) - (α^2 + 2 * ε * λ * α) / (_ξ_^2)))), ξ, Inf, rtol=1e-5)[1]




#Definition of U_λ 
U_λ(ξ, λ) = (1 - 2 / ξ) * (1 + λ^2 / ξ^2)

#Defiition of  R  
function R̃(ξ, ε, α, λ, ϵ_σ)
    temp = ε^2 - U_λ(ξ, λ) - α * (α + 2 * ε * λ * ϵ_σ) / ξ^2
    if temp < 0
        return 0
    else
        return temp
    end
    #println("R = ", temp, " for ε = ",ε, " and λ = ", λ, " ksi = ", ksi )
    #return temp

end

#Definicja S w https://www.actaphys.uj.edu.pl/fulltext?series=Sup&vol=16&aid=6-A13 poniżej wzorów (11). 
S(ξ, ε, λ, α, ϵ_σ, ϵ_r) = exp(-(ε + ϵ_σ * v * sqrt(ε^2 - 1) * sin(φ - ϵ_σ * ϵ_r * (π / 2 - X(ξ, ε, λ, α)))) * β / sqrt(1 - v^2))

function λ_max(ξ, ε, α, ϵ_σ)
    if ξ == 2 && sign(α) * ϵ_σ == 1
        temp = -(α^2 - 4 * ε^2) / (2 * ε * α)
    else
        temp = ξ / (ξ - 2) * (-ϵ_σ * ε * α + sqrt(α^2 * (ε^2 + 2 / ξ - 1) + (ξ - 2) * (ξ * (ε^2 - 1) + 2)))
    end
    return temp
end

function λ_c(α, ε, ϵ_σ, ξ)
    # println("alfa = ", α, " ε =  ", ε, "ϵ_σ = " ,ϵ_σ, " ξ = ", ξ)
    poly_coeffs = [-α^4 - α^6 * (-1 + ε^2),
        (-4 * ϵ_σ * α^3 * ε - 6 * ϵ_σ * α^5 * ε * (-1 + ε^2)),
        (16 - 2 * α^2 - 4 * ϵ_σ^2 * α^2 * ε^2 + 18 * α^2 * (-1 + ε^2) - 3 * α^4 * (-1 + ε^2) - 12 * ϵ_σ^2 * α^4 * ε^2 * (-1 + ε^2)),
        (-4 * ϵ_σ * α * ε + 36 * ϵ_σ * α * ε * (-1 + ε^2) - 12 * ϵ_σ * α^3 * ε * (-1 + ε^2) - 8 * ϵ_σ^3 * α^3 * ε^3 * (-1 + ε^2)),
        (-1 + 18 * (-1 + ε^2) - 3 * α^2 * (-1 + ε^2) - 12 * ϵ_σ^2 * α^2 * ε^2 * (-1 + ε^2) + 27 * (-1 + ε^2)^2),
        6 * ϵ_σ * α * ε * (-1 + ε^2),
        1 - ε^2]
    poly = Polynomial(poly_coeffs)
    #println("polynomial = ", poly)
    sols_imaginary = PolynomialRoots.roots(Float64.(poly_coeffs))
    #print("sols, imag = ", sols_imaginary)
    sols = filter(sol -> abs(imag(sol)) < 1e-15, sols_imaginary)
    limit_λ = -ϵ_σ * α + 2 + 2 * sqrt(1 - ϵ_σ * α)
    #println("sols = ", sols)
    #println("granica = ", limit_λ)
    if isempty(sols)
        # println("empty")
        return nothing
    else
        #println("else")
        real_sols = broadcast(abs, (real.(sols)))
        maks = maximum(real_sols)
        #print("lambdac = ", maks)
        return maks >= limit_λ ? maks : nothing
    end
end


#functions which are in integrals of J currents
function __jt_integrals__(ξ, λ, ε, α, ϵ_σ, ϵ_r=-1)
    if R̃(ξ, ε, α, λ, ϵ_σ) == 0
        return 0
    else
        # temp = ε*S(ξ, ε, λ, α, ϵ_σ, ϵ_r) /sqrt(R̃(ξ, ε,α, λ,ϵ_σ))  
        temp = ε * S(ξ, ε, λ, α, ϵ_σ, ϵ_r) / sqrt(R̃(ξ, ε, α, λ, ϵ_σ)) / (M * ξ) # Tu powinna być suma 2 składników z  ϵ_σ=+1 i ϵ_σ=-1
        return temp
    end
end


#println("Jt_integrand=",__jt_integrals__(8, 1/2, 2, 1/111, 1, -1) )


function __jr_integrals__(ξ, λ, ε, α, ϵ_σ, ϵ_r=-1)
    if R̃(ξ, ε, α, λ, ϵ_σ) == 0
        return 0
    else
        temp = ϵ_r * S(ξ, ε, λ, α, ϵ_σ, ϵ_r)
        return temp
    end
end
function __jφ_integrals__(ξ, λ, ε, α, ϵ_σ, ϵ_r=-1)
    #println("test")
    if R̃(ξ, ε, α, λ, ϵ_σ) == 0
        return 0
    else
        temp = ϵ_σ * S(ξ, ε, λ, α, ϵ_σ, ϵ_r) * (λ + ϵ_σ * α * ε) / sqrt(R̃(ξ, ε, α, λ, ϵ_σ))
        #println("result = ", temp)
        return temp
    end
end

#Calculation of integrals in different J current components ABS
function jt_ABS_integrals(f, ksi, alfa, m_0)
    temp1(λ, ε) = -m_0^3 / ksi * f(ksi, λ, ε, alfa, 1, -1) #eps_sigma = 1; eps_r = -1
    temp2(λ, ε) = -m_0^3 / ksi * f(ksi, λ, ε, alfa, -1, -1) #eps_sigma = -1; eps_r = -1
    result1, err1 = quadgk(ε -> quadgk(λ -> temp1(λ, ε), 0, λ_c(alfa, ε, 1, ksi))[1], 1, Inf)
    result2, err2 = quadgk(ε -> quadgk(λ -> temp1(λ, ε), 0, λ_c(alfa, ε, -1, ksi))[1], 1, Inf)
    result = result1 + result2
    return result

end

function jφ_ABS_integrals(f, ksi, alfa, m_0, M)
    temp1(λ, ε) = M * m_0^3 / ksi * f(ksi, λ, ε, alfa, 1, -1) #eps_sigma = 1; eps_r = -1
    temp2(λ, ε) = M * m_0^3 / ksi * f(ksi, λ, ε, alfa, -1, -1) #eps_sigma = -1; eps_r = -1
    result1, err1 = quadgk(ε -> quadgk(λ -> temp1(λ, ε), 0, λ_c(alfa, ε, 1, ksi))[1], 1, Inf)
    result2, err2 = quadgk(ε -> quadgk(λ -> temp2(λ, ε), 0, λ_c(alfa, ε, -1, ksi))[1], 1, Inf)
    result = result1 + result2
    return result
end
function jr_ABS_integrals(f, ksi, alfa, m_0)
    temp1(λ, ε) = (m_0^3 * ksi) / (ksi * (ksi - 2) + alfa^2) * (f(ksi, λ, ε, alfa, 1, -1)) #eps_sigma = 1; eps_r = -1
    temp2(λ, ε) = (m_0^3 * ksi) / (ksi * (ksi - 2) + alfa^2) * (f(ksi, λ, ε, alfa, -1, -1)) #eps_sigma = -1; eps_r = -1
    result1, err = quadgk(ε -> quadgk(λ -> temp1(λ, ε), 0, λ_c(alfa, ε, 1, ksi))[1], 1, Inf)
    result2, err = quadgk(ε -> quadgk(λ -> temp2(λ, ε), 0, λ_c(alfa, ε, -1, ksi))[1], 1, Inf)
    result = result1 + result2
    return result
end
function ξ_ph(eps_sigma, alfa)
    temp = 2 + 2 * cos(2 / 3 * acos(-eps_sigma * alfa))
    # println("ξ_ph = ", temp)
    return temp
end
function ξ_mb(eps_sigma, alfa)

    temp = 2 - eps_sigma * alfa + 2 * sqrt(1 - eps_sigma * alfa)
    #println("ξ_mb = ", temp)
    return temp
end

function ε_min(ξ, α, ϵ_σ)
    if ξ < ξ_ph(ϵ_σ, α)
        temp = Inf
    elseif ξ_ph(ϵ_σ, α) < ξ && ξ < ξ_mb(ϵ_σ, α)
        temp = sqrt((-2 * ϵ_σ * α * (α^2 + (ξ - 2) * ξ) * ξ^(-1 / 2) + α^2 * (5 - 3 * ξ) + (ξ - 3) * (ξ - 2)^2 * ξ) / (ξ * ((ξ - 3)^2 * ξ - 4 * α^2)))
    elseif ξ >= ξ_mb(ϵ_σ, α)
        temp = 1
    end
    # print("eps_min = ", temp)
    return temp
end
#Calculation of integrals in different J current components SCATT
function jt_SCATT_integrals(f, ksi, alfa, m_0)
    temp1(λ, ε) = -m_0^3 / ksi * (f(ksi, λ, ε, alfa, 1, 1)) #eps_sigma = 1; eps_r = 1
    temp2(λ, ε) = -m_0^3 / ksi * (f(ksi, λ, ε, alfa, -1, 1)) #eps_sigma = -1; eps_r = 1
    temp3(λ, ε) = -m_0^3 / ksi * (f(ksi, λ, ε, alfa, 1, -1)) #eps_sigma = 1; eps_r = -1
    temp4(λ, ε) = -m_0^3 / ksi * (f(ksi, λ, ε, alfa, -1, -1)) #eps_sigma = 1; eps_r = 1
    result1, err1 = quadgk(ε -> quadgk(λ -> temp1(λ, ε), λ_c(alfa, ε, 1, ksi), λ_max(ksi, ε, alfa, 1))[1], ε_min(ksi, alfa, 1), Inf) #lower boundary = λ_c; upper_boundary = λ_max 
    result2, err2 = quadgk(ε -> quadgk(λ -> temp2(λ, ε), λ_c(alfa, ε, -1, ksi), λ_max(ksi, ε, alfa, -1))[1], ε_min(ksi, alfa, -1), Inf)
    result3, err3 = quadgk(ε -> quadgk(λ -> temp3(λ, ε), λ_c(alfa, ε, 1, ksi), λ_max(ksi, ε, alfa, 1))[1], ε_min(ksi, alfa, 1), Inf)
    result4, err4 = quadgk(ε -> quadgk(λ -> temp4(λ, ε), λ_c(alfa, ε, -1, ksi), λ_max(ksi, ε, alfa, -1))[1], ε_min(ksi, alfa, -1), Inf)
    result = result1 + result2 + result3 + result4
    return result
end
function jφ_SCATT_integrals(f, ksi, alfa, m_0, M)
    temp1(λ, ε) = M * m_0^3 / ksi * (f(ksi, λ, ε, alfa, 1, 1)) #eps_sigma = 1; eps_r = 1
    temp2(λ, ε) = M * m_0^3 / ksi * (f(ksi, λ, ε, alfa, -1, 1)) #eps_sigma = -1; eps_r = 1
    temp3(λ, ε) = M * m_0^3 / ksi * (f(ksi, λ, ε, alfa, 1, -1)) #eps_sigma = 1; eps_r = -1
    temp4(λ, ε) = M * m_0^3 / ksi * (f(ksi, λ, ε, alfa, -1, -1)) #eps_sigma = 1; eps_r = 1t
    result1, err1 = quadgk(ε -> quadgk(λ -> temp1(λ, ε), λ_c(alfa, ε, 1, ksi), λ_max(ksi, ε, alfa, 1))[1], ε_min(ksi, alfa, 1), Inf) #lower boundary = λ_c; upper_boundary = λ_max 
    result2, err2 = quadgk(ε -> quadgk(λ -> temp2(λ, ε), λ_c(alfa, ε, -1, ksi), λ_max(ksi, ε, alfa, -1))[1], ε_min(ksi, alfa, -1), Inf)
    result3, err3 = quadgk(ε -> quadgk(λ -> temp3(λ, ε), λ_c(alfa, ε, 1, ksi), λ_max(ksi, ε, alfa, 1))[1], ε_min(ksi, alfa, 1), Inf)
    result4, err4 = quadgk(ε -> quadgk(λ -> temp4(λ, ε), λ_c(alfa, ε, -1, ksi), λ_max(ksi, ε, alfa, -1))[1], ε_min(ksi, alfa, -1), Inf)
    result = result1 + result2 + result3 + result4
    return result
end
function jr_SCATT_integrals(f, ksi, alfa, m_0)
    temp1(λ, ε) = m_0^3 * ksi / (ksi * (ksi - 2) + alfa^2) * (f(ksi, λ, ε, alfa, 1, 1)) #eps_sigma = 1; eps_r = 1
    temp2(λ, ε) = m_0^3 * ksi / (ksi * (ksi - 2) + alfa^2) * (f(ksi, λ, ε, alfa, -1, 1)) #eps_sigma = -1; eps_r = 1
    temp3(λ, ε) = m_0^3 * ksi / (ksi * (ksi - 2) + alfa^2) * (f(ksi, λ, ε, alfa, 1, -1)) #eps_sigma = 1; eps_r = -1
    temp4(λ, ε) = m_0^3 * ksi / (ksi * (ksi - 2) + alfa^2) * (f(ksi, λ, ε, alfa, -1, -1)) #eps_sigma = 1; eps_r = 1
    result1, err1 = quadgk(ε -> quadgk(λ -> temp1(λ, ε), λ_c(alfa, ε, 1, ksi), λ_max(ksi, ε, alfa, 1))[1], ε_min(ksi, alfa, 1), Inf) #lower boundary = λ_c; upper_boundary = λ_max 
    result2, err2 = quadgk(ε -> quadgk(λ -> temp2(λ, ε), λ_c(alfa, ε, -1, ksi), λ_max(ksi, ε, alfa, -1))[1], ε_min(ksi, alfa, -1), Inf)
    result3, err3 = quadgk(ε -> quadgk(λ -> temp3(λ, ε), λ_c(alfa, ε, 1, ksi), λ_max(ksi, ε, alfa, 1))[1], ε_min(ksi, alfa, 1), Inf)
    result4, err4 = quadgk(ε -> quadgk(λ -> temp4(λ, ε), λ_c(alfa, ε, -1, ksi), λ_max(ksi, ε, alfa, -1))[1], ε_min(ksi, alfa, -1), Inf)
    result = result1 + result2 + result3 + result4
    return result
end
function create_r_tbl(start, ending, step, M)
    seq1 = collect(start:step:ending)
    seq2 = [2^k for k in 1:10]
    seq = unique(vcat(seq1, seq2)) #unique - removes duplicates; vcat - adds sequances vertically
    result = seq .* M#filter(x -> !(x in [xi_hor, xi_ph, xi_mb]), seq) # .* operator multiplies sequence by a number elementwise
    return result
end
############AFTER DEFINITIONS################
ksi = 10#abs(rand(3:20))
M = 1;
m_0 = 1; #M - mass of the black hole, m_0 - mass of the particles
φ_table = LinRange(-π, π, 73)
r_table = create_r_tbl(2, 18, 1 / 16, M)
φ = 1 #some "random" angle
for i in 1:1

    alfa = 0.0001
    println("ksi = ", ksi)
    # try
    J_t_ABS = jt_ABS_integrals(__jt_integrals__,ksi,alfa,m_0)
    J_r_ABS = jr_ABS_integrals(__jr_integrals__, ksi, alfa, m_0)
    J_φ_ABS = jφ_ABS_integrals(__jφ_integrals__, ksi, alfa, m_0, M)
    ###################################################
    # J_t_SCATT = jt_SCATT_integrals(__jt_integrals__, ksi, alfa, m_0)
    # J_r_SCATT= jr_SCATT_integrals(__jr_integrals__, ksi, alfa, m_0)
    # J_φ_SCATT = jφ_SCATT_integrals(__jφ_integrals__, ksi, alfa, m_0,M)

    println("J_t_ABS = ", J_t_ABS)
    println("J_r_ABS = ", J_r_ABS)
    println("J_φ_ABS = ", J_φ_ABS)
    #
    #println("J_t_SCATT = ", J_t_SCATT)
    #println("J_r_SCATT = ", J_r_SCATT)
    #println("J_φ_SCATT = ", J_φ_SCATT)
    ########################################################## UWAGA! Tutaj ksi powinno byc r-em, zmienic!!!! W kodzie prof. Odrzywołka jest to w jednostkach masy!!!!
    J_X_ABS = J_r_ABS * cos(φ) - J_φ_ABS * ksi * M * sin(φ) #J^x  = J^r Cosφ - (J^φ) r Sinφ 
    J_Y_ABS = J_r_ABS * sin(φ) - J_φ_ABS * ksi * M * cos(φ) #J^y  = J^r Sinφ + J^φ r Cosφ
    println("J_x_ABS = ", J_X_ABS)
    println("J_y_ABS = ", J_Y_ABS)

    #catch e_rror
    #   println(e_rror)
    #end
end
