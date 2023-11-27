#= 
Author: Ksymena Poradzisz
Contact: ksymena.poradzisz@gmail.com
Affiliation: Jagiellonian University
Created: [2023-10-21]
Updated: [2023-11-18]
Description:
This Julia script is intended to compute integrals for black hole accretion currents J^\mu
=#
#using HCubature
using Symbolics
using QuadGK
using DoubleExponentialFormulas, CSV
@variables M, m_0, ξ, ϵ_σ, ϵ_r, ε, α, λ, _ξ_,φ
machine_epsilon = eps() #the smallest number which can be added to 1.0
#X(ξ, ε, λ, α) = quadde(_ξ_-> (λ + (α*ε)/(1-2/_ξ_))/( ( (α^2)/(1-2/_ξ_) + _ξ_^2  ) * sqrt(ε^2 - (1-2/_ξ_)*(1+λ^2/_ξ_^2) - (α^2+2*ε*λ*α)/(_ξ_^2)) ), ξ, Inf, rtol = 1e-5)[1] 
X(ξ, ε, λ, α) = quadde(_ξ_-> (λ + (α*ε)/(1-2/_ξ_))/( ( (α^2)/(1-2/_ξ_) + _ξ_^2  ) * sqrt(ε^2 - (1-2/_ξ_)*(1+λ^2/_ξ_^2) - (α^2+2*ε*λ*α)/(_ξ_^2)) ), ξ, Inf, rtol = 1e-5)[1] 

v = 0.5 # predkosc w jedn. c
β = 2   # Thermodynamic β=1/(kT), aka coolness
φ = 1   # some "random" polar angle

λ_c = sqrt(12/(1- (4)/(3*ε)/sqrt(9*ε^2)+1)) # 


#Definicja U_λ TODO
U_λ(ξ, λ) = (1-2/ξ)*(1+λ^2/ξ^2)

#Definicja R  TODO
R̃(ξ,ε,α,λ ) = ε^2 - U_λ(ξ,λ) - α * (α + 2*ε*λ*ϵ_σ)/ξ^2

#Definicja S w https://www.actaphys.uj.edu.pl/fulltext?series=Sup&vol=16&aid=6-A13 poniżej wzorów (11). 
S(ξ, ε, λ, α, ϵ_σ, ϵ_r)  = exp(-(ε+ϵ_σ*v*sqrt(ε^2-1)*sin( (-X(ξ,ε, λ, α))*(φ-ϵ_σ*ϵ_r)))*β/sqrt(1-v^2) )


#S(ξ, ε, λ, α, ϵ_σ)  = exp(-(ε+ϵ_σ*v*sqrt(ε^2-1)*sin( (-sin(λ))*(φ-ϵ_σ*ϵ_r)))*(β)/sqrt(1-v^2) )
#S(ξ, ε, λ, α, ϵ_σ)  = exp(sin( (-sin(λ))*(φ-ϵ_σ))*β/sqrt(1-v^2))

#_j_t_integrals_(ξελαϵ_σ::SVector(5, Float64))  =  ξελαϵ_σ[2]* S(ξελαϵ_σ[1],ξελαϵ_σ[2],ξελαϵ_σ[3],ξελαϵ_σ[4],ξελαϵ_σ[5]) / sqrt( R̃(ξελαϵ_σ[1],ξελαϵ_σ[2],ξελαϵ_σ[4],ξελαϵ_σ[3])) 
function __jt_integrals__(ξ, λ, ε, α, ϵ_σ, ϵ_r  = -1)
    temp = ε*S(ξ, ε, λ, α, ϵ_σ, ϵ_r) /sqrt(R̃(ξ, ε,α, λ))  
    return temp  
    
end

function jt_integrals(f, ksi, alfa, eps_sigma, eps_r)
    temp(λ, ε) = f(ksi, λ, ε, alfa, eps_sigma, eps_r)
    result, err = quadgk(ε->quadgk(λ->temp(λ, ε), 0, sqrt(12/(1- (4)/(3*ε)/sqrt(9*ε^2)+1)))[1], 1, Inf)
    return result
    
end


function __jr_integrals__(ξ,  λ,ε, α, ϵ_σ, ϵ_r=-1)
    temp = ϵ_r * S(ξ, ε, λ, α, ϵ_σ, ϵ_r) 
    return temp
end


function jr_integrals(f,  ksi, alfa, eps_sigma, eps_r)
    temp(λ, ε) = f(ksi,  λ,ε, alfa, eps_sigma, eps_r)
    result, err = quadgk( ε->quadgk(λ-> temp(λ, ε), 0, sqrt(12/(1- (4)/(3*ε)/sqrt(9*ε^2)+1)))[1], 1, Inf)
   # result, err = hcubature(s->hcubature(r-> temp(r[1], s[1]), [0], [sqrt(12/(1- (4)/(3*s[1])/sqrt(9*s[1]^2)+1))])[1],[1],[10.0^18])
    return result
end


ksi = abs(rand(9:20)) 

for i in 1:1
    
    alfa = 0
    eps_sigma = 1
    eps_r = -1
    try
        res = jr_integrals(__jr_integrals__,  ksi, alfa, eps_sigma, eps_r)
        println(res)
    catch e_rror
        println(e_rror)
    end
end












