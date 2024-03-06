#= 
Author: Ksymena Poradzisz
Contact: ksymena.poradzisz@gmail.com
Affiliation: Jagiellonian University
Created: [2023-10-21]
Updated: [2024-03-06]
Description:
This Julia script is intended to compute integrals for black hole accretion currents J^\mu
=#
using Symbolics
using QuadGK
using Roots
#using DoubleExponentialFormulas, CSV
@variables M, m_0, ξ, ϵ_σ, ϵ_r, ε, α, λ, _ξ_,φ
machine_epsilon = eps() #the smallest number which can be added to 1.0
v = 0.5 # predkosc w jedn. c
β = 2   # Thermodynamic β=1/(kT), aka coolness
φ = 1   # some "random" polar angle

################FUNCTIONS############

X(ξ, ε, λ, α) = quadgk(_ξ_-> (λ + (α*ε)/(1-2/_ξ_))/( ( (α^2)/(1-2/_ξ_) + _ξ_^2  ) * sqrt(Complex(ε^2 - (1-2/_ξ_)*(1+λ^2/_ξ_^2) - (α^2+2*ε*λ*α)/(_ξ_^2))) ), ξ, Inf, rtol = 1e-5)[1] 

#Definicja U_λ 
U_λ(ξ, λ) = (1-2/ξ)*(1+λ^2/ξ^2)

#Definicja R  
function R̃(ξ,ε,α,λ,ϵ_σ )
	temp = ε^2 - U_λ(ξ,λ) - α * (α + 2*ε*λ*ϵ_σ)/ξ^2
	return temp

end

#Definicja S w https://www.actaphys.uj.edu.pl/fulltext?series=Sup&vol=16&aid=6-A13 poniżej wzorów (11). 
S(ξ, ε, λ, α, ϵ_σ, ϵ_r)  = exp(-(ε+ϵ_σ*v*sqrt(ε^2-1)*sin( (π/2 -X(ξ,ε, λ, α))*(φ-ϵ_σ*ϵ_r)))*β/sqrt(1-v^2) )

function λ_max(ξ, ε, α,ϵ_σ )
    if ξ == 2 && sign(α)*ϵ_σ == 1
        temp = -(α^2-4*ε^2)/(2*ε*α)
    else
        temp = ξ/(ξ-2) * (-ϵ_σ*ε*α + sqrt(α^2*(ε^2+2/ξ - 1)+(ξ-2)*(ξ*(ε^2-1)+2)))
    end
    return temp
end

function λ_c(α,ε, ϵ_σ, ξ)
    function eq_func(λ)
        Y_temp = α^2 + 2*ε*α*ϵ_σ*λ
        eq = (ε^2 - 1) * (27 * (ε^2 - 1) * λ^4 - Y_temp^3 + 18 * λ^2 * Y_temp) + (16 * λ^2 - Y_temp^2)
        return eq
    end
    limit_λ = -ε*α + 2 + 2*sqrt(1-ε*α)
    sols = find_zero(eq_func, (limit_λ, λ_max(ξ,ε, α, ϵ_σ )))
    print(sols)
    maks = maximum(sol for sol in sols if sol >= limit_λ)
    return maks
end


#functions which are in integrals of J currents
function __jt_integrals__(ξ, λ, ε, α, ϵ_σ, ϵ_r  = -1)
    if R̃(ξ, ε,α, λ,ϵ_σ) < 0
       # println("R̃(ξ, ε,α, λ,ϵ_σ) = ", R̃(ξ, ε,α, λ,ϵ_σ), "for ξ = ", ξ, " ε = ", ε, " λ = ", λ)
       # println("Lambda max = ", λ_max(ε,ξ))
    end
    temp = ε*S(ξ, ε, λ, α, ϵ_σ, ϵ_r) /sqrt(R̃(ξ, ε,α, λ,ϵ_σ))  
    return temp  
    
end
function __jr_integrals__(ξ,  λ,ε, α, ϵ_σ, ϵ_r=-1)
    temp = ϵ_r * S(ξ, ε, λ, α, ϵ_σ, ϵ_r) 
    return temp
end
function __jφ_integrals__(ξ, λ,ε, α, ϵ_σ, ϵ_r = -1 )
    temp = ϵ_σ * S(ξ, ε, λ, α, ϵ_σ, ϵ_r)*(λ + ϵ_σ * α * ε)/sqrt(R̃(ξ, ε, α, λ,ϵ_σ))
    return temp
end

#Calculation of integrals in different J current components ABS
function jt_ABS_integrals(f, ksi, alfa, eps_sigma, eps_r)
    temp(λ, ε) = f(ksi, λ, ε, alfa, eps_sigma, eps_r)
    result, err = quadgk(ε->quadgk(λ->temp(λ, ε), 0, sqrt(12/(1- (4)/(3*ε)/sqrt(9*ε^2)+1)))[1], 1, Inf)
    return result
    
end

function jφ_ABS_integrals(f, ksi, alfa, eps_sigma, eps_r)
    temp(λ,ε) = f(ksi, λ, ε, alfa, eps_sigma, eps_r)
    result,err = quadgk(ε-> quadgk(λ-> temp(λ, ε), 0, sqrt(12/(1- (4)/(3*ε)/sqrt(9*ε^2)+1)))[1], 1, Inf)
    return result
end
function jr_ABS_integrals(f,  ksi, alfa, eps_sigma, eps_r)
    temp(λ, ε) = f(ksi,  λ,ε, alfa, eps_sigma, eps_r)
    result, err = quadgk( ε->quadgk(λ-> temp(λ, ε), 0, sqrt(12/(1- (4)/(3*ε)/sqrt(9*ε^2)+1)))[1], 1, Inf)
    return result
end
function ε_min(ξ)
    if ξ<=3
        return Inf
    elseif ξ > 3 && ξ < 4
        return sqrt((1-2/ξ)*(1+1/(ξ-3)))
    elseif ξ >=4
        return 1
    end
end
#Calculation of integrals in different J current components SCATT
function jt_SCATT_integrals(f, ksi, alfa, eps_sigma, eps_r)
    temp(λ, ε) = f(ksi,  λ,ε, alfa, eps_sigma, eps_r)
    #λ_c(ε) = sqrt(12/(1- (4)/(3*ε)/sqrt(9*ε^2)+1))
    #λ_max(ε,ξ) = ξ *sqrt( (ε^2)/(1-2/ξ) -1)   
    result, err = quadgk(ε-> quadgk(λ->temp(λ, ε), λ_c(ε), λ_max(ε, ksi))[1] , 1, Inf) #lower boundary = λ_c; upper_boundary = λ_max 
    return result
end
function jφ_SCATT_integrals(f, ksi, alfa, eps_sigma, eps_r)
    temp(λ,ε) = f(ksi, λ,ε, alfa, eps_sigma, eps_r)
    #λ_c(ε) = sqrt(12/(1- (4)/(3*ε)/sqrt(9*ε^2)+1))
    #λ_max(ε,ξ) = ξ *sqrt( (ε^2)/(1-2/ξ) -1) 
    result, err = quadgk(ε-> quadgk(λ->temp(λ, ε), λ_c(ε), λ_max(ε, ksi))[1], ε_min(ksi), Inf)
    return result
end
function jr_SCATT_integrals(f, ksi, alfa, eps_sigma, eps_r)
    temp(λ,ε) = f(ksi, λ,ε, alfa, eps_sigma, eps_r)
   # λ_c(ε) = sqrt(12/(1- (4)/(3*ε)/sqrt(9*ε^2)+1))
    #λ_max(ε,ξ) = ξ *sqrt( (ε^2)/(1-2/ξ) -1) 
    result, err = quadgk(ε-> quadgk(λ->temp(λ, ε), λ_c(ε), λ_max(ε, ksi))[1], ε_min(ksi), Inf)
    return result
end
ksi = abs(rand(10:20)) 

for i in 1:1
    
    alfa = 0.0001
    eps_sigma = 1
    eps_r = -1
    print(λ_c(alfa, 2, eps_sigma, ksi))
    #test_phi = jφ_SCATT_integrals(__jφ_integrals__, ksi, alfa, eps_sigma, eps_r)
    #println("Phi SCATT = ",test_phi)
    try
        #test_t = jt_SCATT_integrals(__jt_integrals__, ksi, alfa, eps_sigma, eps_r)
       # println("Time SCATT = ", test_t)
    catch e_rror
        println(e_rror)
    end
end
