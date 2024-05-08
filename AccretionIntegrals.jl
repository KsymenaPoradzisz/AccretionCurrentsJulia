#= 
Author: Ksymena Poradzisz
Contact: ksymena.poradzisz@gmail.com
Affiliation: Jagiellonian University
Created: [2023-10-21]
Updated: [2024-05-07]
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
3. Run code from Linux shell: julia AccretionIntegrals_20240507.jl 

Expected outcome:

 

=#

using Symbolics
using QuadGK
using PolynomialRoots
using Polynomials
#using DoubleExponentialFormulas, CSV
@variables M, m_0, ξ, ϵ_σ, ϵ_r, ε, α, λ, _ξ_,φ
machine_epsilon = eps() #the smallest number which can be added to 1.0
v = 0.5 # predkosc w jedn. c
β = 2   # Thermodynamic β=1/(kT), aka coolness
φ = 1   # some "random" polar angle
################FUNCTIONS############

X(ξ, ε, λ, α) = quadgk(_ξ_-> (λ + (α*ε)/(1-2/_ξ_))/( ( (α^2)/(1-2/_ξ_) + _ξ_^2  ) * sqrt(Complex(ε^2 - (1-2/_ξ_)*(1+λ^2/_ξ_^2) - (α^2+2*ε*λ*α)/(_ξ_^2))) ), ξ, Inf, rtol = 1e-5)[1] 


println("X=",X(11,2,1/2,1/111))

#Definicja U_λ 
U_λ(ξ, λ) = (1-2/ξ)*(1+λ^2/ξ^2)

#Definicja R  
function R̃(ξ,ε,α,λ,ϵ_σ )
	temp = ε^2 - U_λ(ξ,λ) - α * (α + 2*ε*λ*ϵ_σ)/ξ^2
    if temp < 0
        return 0
    else 
        return temp
    end
	#println("R = ", temp, " for ε = ",ε, " and λ = ", λ, " ksi = ", ksi )
    #return temp

end

#Definicja S w https://www.actaphys.uj.edu.pl/fulltext?series=Sup&vol=16&aid=6-A13 poniżej wzorów (11). 
#S(ξ, ε, λ, α, ϵ_σ, ϵ_r)  = exp(-(ε+ϵ_σ*v*sqrt(ε^2-1)*sin( (π/2 -X(ξ,ε, λ, α))*(φ-ϵ_σ*ϵ_r)))*β/sqrt(1-v^2) )
S(ξ, ε, λ, α, ϵ_σ, ϵ_r)  = exp(-(ε+ϵ_σ*v*sqrt(ε^2-1)*sin(φ-ϵ_σ*ϵ_r*(π/2-X(ξ,ε, λ, α))))*β/sqrt(1-v^2) )

function λ_max(ξ, ε, α,ϵ_σ )
    if ξ == 2 && sign(α)*ϵ_σ == 1
        temp = -(α^2-4*ε^2)/(2*ε*α)
    else
        temp = ξ/(ξ-2) * (-ϵ_σ*ε*α + sqrt(α^2*(ε^2+2/ξ - 1)+(ξ-2)*(ξ*(ε^2-1)+2)))
    end
    return temp
end
#Porównując kod z mathematica, wartość λ_c jest liczona poprawnie. 
function λ_c(α, ε, ϵ_σ, ξ)
   # println("alfa = ", α, " ε =  ", ε, "ϵ_σ = " ,ϵ_σ, " ξ = ", ξ)
    poly_coeffs = [-α^4 - α^6 * (-1 + ε^2),
                       (-4 * ϵ_σ * α^3 * ε - 6 * ϵ_σ * α^5 * ε * (-1 + ε^2)),
                       (16 - 2 * α^2 - 4 * ϵ_σ^2 * α^2 * ε^2 + 18 * α^2 * (-1 + ε^2) - 3 * α^4 * (-1 + ε^2) - 12 * ϵ_σ^2 * α^4 * ε^2 * (-1 + ε^2)),
                       (-4 * ϵ_σ * α * ε + 36 * ϵ_σ * α * ε * (-1 + ε^2) - 12 * ϵ_σ * α^3 * ε * (-1 + ε^2) - 8 * ϵ_σ^3 * α^3 * ε^3 * (-1 + ε^2)),
                       (-1 + 18 * (-1 + ε^2) - 3 * α^2 * (-1 + ε^2) -  12 * ϵ_σ^2 * α^2 * ε^2 * (-1 + ε^2) + 27 * (-1 + ε^2)^2),
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
        real_sols = broadcast(abs,(real.(sols)))
        maks = maximum(real_sols)
        #print("lambdac = ", maks)
        return maks >= limit_λ ? maks : nothing
    end
end


#functions which are in integrals of J currents
function __jt_integrals__(ξ, λ, ε, α, ϵ_σ, ϵ_r  = -1)
    if R̃(ξ, ε,α, λ,ϵ_σ) == 0
        return 0
    else
        # temp = ε*S(ξ, ε, λ, α, ϵ_σ, ϵ_r) /sqrt(R̃(ξ, ε,α, λ,ϵ_σ))  
        temp = ε*S(ξ, ε, λ, α, ϵ_σ, ϵ_r) /sqrt(R̃(ξ, ε,α, λ,ϵ_σ))/(M*ξ) # Tu powinna być suma 2 składników z  ϵ_σ=+1 i ϵ_σ=-1
        return temp  
    end 
end


println("Jt_integrand=",__jt_integrals__(8, 1/2, 2, 1/111, 1, -1) )


function __jr_integrals__(ξ,  λ,ε, α, ϵ_σ, ϵ_r=-1)
    if R̃(ξ, ε,α, λ,ϵ_σ) == 0
        return 0
    else
        temp = ϵ_r * S(ξ, ε, λ, α, ϵ_σ, ϵ_r) 
        return temp
    end
end
function __jφ_integrals__(ξ, λ,ε, α, ϵ_σ, ϵ_r = -1 )
    #println("test")
    if R̃(ξ, ε,α, λ,ϵ_σ) == 0
        return 0
    else
        temp = ϵ_σ * S(ξ, ε, λ, α, ϵ_σ, ϵ_r)*(λ + ϵ_σ * α * ε)/sqrt(R̃(ξ, ε, α, λ,ϵ_σ))
        #println("result = ", temp)
        return temp
    end
end

#Calculation of integrals in different J current components ABS
function jt_ABS_integrals(f, ksi, alfa, eps_sigma, eps_r)
    temp(λ, ε) = f(ksi, λ, ε, alfa, eps_sigma, eps_r)
    result, err = quadgk(ε->quadgk(λ->temp(λ, ε), 0, λ_c(alfa, ε, eps_sigma, ksi))[1], 1, Inf)
    return result
    
end

function jφ_ABS_integrals(f, ksi, alfa, eps_sigma, eps_r) #ta funkcja chyba działa
    temp(λ,ε) = f(ksi, λ, ε, alfa, eps_sigma, eps_r)
    result,err = quadgk(ε-> quadgk(λ-> temp(λ, ε), 0, λ_c(alfa, ε, eps_sigma, ksi))[1], 1, Inf)
    return result
end
function jr_ABS_integrals(f,  ksi, alfa, eps_sigma, eps_r)
    temp(λ, ε) = f(ksi,  λ,ε, alfa, eps_sigma, eps_r)
    result, err = quadgk( ε->quadgk(λ-> temp(λ, ε), 0, λ_c(alfa, ε, eps_sigma, ksi))[1], 1, Inf)
    return result
end
function ξ_ph(eps_sigma, alfa)
    temp = 2 + 2*cos(2/3*acos(-eps_sigma*alfa))
   # println("ξ_ph = ", temp)
    return temp
end
function ξ_mb(eps_sigma, alfa)
  
    temp = 2 - eps_sigma*alfa + 2*sqrt(1-eps_sigma*alfa)
    #println("ξ_mb = ", temp)
    return temp 
end

function ε_min(ξ, α, ϵ_σ)
    if ξ < ξ_ph(ϵ_σ, α)
        temp = Inf
    elseif ξ_ph(ϵ_σ, α) < ξ && ξ < ξ_mb(ϵ_σ, α)
        temp = sqrt((-2 * ϵ_σ * α * (α^2 + (ξ - 2) * ξ) * ξ^(-1/2) + α^2 * (5 - 3*ξ) + (ξ - 3) * (ξ - 2)^2 * ξ) / (ξ * ((ξ - 3)^2 * ξ - 4*α^2)))
    elseif ξ >= ξ_mb(ϵ_σ, α)
        temp = 1
    end
    print("eps_min = ", temp)
    return temp
end
#Calculation of integrals in different J current components SCATT
function jt_SCATT_integrals(f, ksi, alfa, eps_sigma, eps_r)
    temp(λ, ε) = f(ksi,  λ,ε, alfa, eps_sigma, eps_r)
    result, err = quadgk(ε-> quadgk(λ->temp(λ, ε), λ_c(alfa, ε, eps_sigma, ksi), λ_max(ksi, ε, alfa, eps_sigma))[1] , ε_min(ksi, alfa, eps_sigma), Inf) #lower boundary = λ_c; upper_boundary = λ_max 
    return result
end
function jφ_SCATT_integrals(f, ksi, alfa, eps_sigma, eps_r)
    temp(λ,ε) = f(ksi, λ,ε, alfa, eps_sigma, eps_r)
    result, err = quadgk(ε-> quadgk(λ->temp(λ, ε), λ_c(alfa, ε, eps_sigma, ksi), λ_max(ksi, ε, alfa, eps_sigma))[1], ε_min(ksi, alfa, eps_sigma), Inf)
    #println("result = ", result)
    return result
end
function jr_SCATT_integrals(f, ksi, alfa, eps_sigma, eps_r)
    temp(λ,ε) = f(ksi, λ,ε, alfa, eps_sigma, eps_r)
    result, err = quadgk(ε-> quadgk(λ->temp(λ, ε), λ_c(alfa, ε, eps_sigma, ksi), λ_max(ksi, ε, alfa, eps_sigma))[1], ε_min(ksi, alfa, eps_sigma), Inf)
    return result
end
ksi = abs(rand(2:20)) 

for i in 1:1
    
    alfa = 0.0001
    eps_sigma = 1
    eps_r = -1
    println("ksi = ", ksi)
    try
    #print("lambda_c = ", λ_c(Float64(alfa), 2., Float64(eps_sigma), Float64(ksi)))
    test_phi_abs = jφ_ABS_integrals(__jφ_integrals__, ksi, alfa, eps_sigma, eps_r)
    test_t_abs =  jt_ABS_integrals(__jt_integrals__, ksi, alfa, eps_sigma, eps_r)
    test_r_abs = jr_ABS_integrals(__jr_integrals__,  ksi, alfa, eps_sigma, eps_r)
    ###########################################################################
    test_phi_scatt  = jφ_SCATT_integrals(__jφ_integrals__, ksi, alfa, eps_sigma, eps_r)
    test_t_scatt = jt_SCATT_integrals(__jt_integrals__, ksi, alfa, eps_sigma, eps_r)
    test_r_scatt = jr_SCATT_integrals(__jr_integrals__, ksi, alfa, eps_sigma, eps_r)
    ########################################################################
    println("Phi ABS = ",test_phi_abs )
    println("t ABS = ", test_t_abs)
    println("r ABS = ", test_r_abs)
    ##############################################
    println("phi SCATT = ", test_phi_scatt)
    println("t SCATT = ", test_t_scatt)
    println("r SCATT = ", test_r_scatt)
    
        #test_t = jt_SCATT_integrals(__jt_integrals__, ksi, alfa, eps_sigma, eps_r)
       # println("Time SCATT = ", test_t)
    catch e_rror
        println(e_rror)
    end
end
