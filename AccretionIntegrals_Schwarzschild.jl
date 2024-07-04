#= 
Author: Ksymena Poradzisz
Contact: ksymena.poradzisz@gmail.com
Affiliation: Jagiellonian University
Created: [2024-06-28]
Updated: [2024-07-04]
Description:
This Julia script is intended to compute integrals for Schwarzschild black hole accretion currents J^\mu 
=#

#using Symbolics
using QuadGK
using DoubleExponentialFormulas
const v = -0.5 # predkosc w jedn. c
const β = 2  # Thermodynamic β=1/(kT), aka coolness
const γ = 1 / sqrt(1 - v^2) #global usage unnecessary
const infinity = 123
M = 1; m_0 = 1;
α_sch = 1 #it is NOT a Kerr parameter - it is related to surface density. Formula (41) in https://arxiv.org/abs/2406.04471
λ_max(ξ, ε) = ξ * sqrt(ε^2 / (1 - 2 / ξ) - 1)

λ_c(ε) = sqrt(12 / (1 - 4 / (1 + 3 * ε / sqrt(9 * ε^2 - 8))^2))
U_λ(ξ, λ) = (1 - 2 / ξ) * (1 + λ^2 / ξ^2)

function R(ξ, ε, λ)
    temp = ε^2 - U_λ(ξ, λ)
    if temp <= 0
        return Inf
    else
        return temp
    end
    #println("R = ", temp, " for ε = ",ε, " and λ = ", λ, " ksi = ", ksi )
    #return temp

end

#X(ξ, ε, λ) = quadgk(_ξ_-> λ/(_ξ_^2*sqrt(R(_ξ_, ε, λ))),ξ,Inf)[1]
X(ξ, ε, λ) = quadde(_ξ_-> λ/(_ξ_^2*sqrt(R(_ξ_, ε, λ))),ξ,Inf)[1]


function ε_min(ξ)
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

#=
function J_t_SCATT(ε, ξ, φ)
    integral= 
               quadgk( λ ->  1*cosh(β*γ*v*sqrt(ε^2-1)*cos(φ)*cos(X(ξ, ε, λ))) *cosh(β*γ*v*sqrt(ε^2-1)*sin(φ)*sin(X(ξ, ε, λ)))/sqrt(R(ξ, ε, λ)), λ_c(ε), λ_max(ξ,ε))[1]
    return integral
end
=#
function J_t_ABS(ξ,φ)
    #coeff = -2*α_sch*m_0^3/ξ
    coeff = 1
    integral = coeff *  quadgk( ε-> ε* quadgk(λ-> cosh(β*γ*v*sqrt(ε^2-1)*sin(φ)*sin(X(ξ, ε, λ)))*
    exp(-β*γ*(ε+v*sqrt(ε^2-1)*cos(φ)*cos(X(ξ,ε,λ))))/sqrt(R(ξ, ε, λ)),0,λ_c(ε))[1],1,infinity)[1]
    return integral
end
function J_t_SCATT(ξ, φ)
    #coeff = -4*α_sch*m_0^3/ξ
    coeff = 1
    integral=  coeff*quadgk( 
                ε-> exp(-β*γ*ε) * ε *
                 quadgk( λ ->  1*cosh(β*γ*v*sqrt(ε^2-1)*cos(φ)*cos(X(ξ, ε, λ))) *cosh(β*γ*v*sqrt(ε^2-1)*sin(φ)*sin(X(ξ, ε, λ)))/sqrt(R(ξ, ε, λ)), λ_c(ε), λ_max(ξ,ε))[1],ε_min(ξ), infinity)[1]
    return integral
end
function J_r_ABS(ξ,φ)
    #coeff = -2*α_sch*m_0^3/(ξ-2) 
    coeff = 1
    integral = coeff *quadgk( ε->  
    quadgk(λ-> cosh(β*γ*v*sqrt(ε^2-1)*sin(φ)*sin(X(ξ, ε, λ)))*exp(-β*γ*(ε+v*sqrt(ε^2-1)*cos(φ)*cos(X(ξ,ε,λ)))),0,λ_c(ε))[1],1,infinity)[1]
    return integral
end
function J_r_SCATT(ξ,φ)
   # coeff = -4*α_sch*m_0^3/(ξ-2)
    coeff = 1
    integral = coeff * quadgk(ε-> exp(-β*γ*ε) * 
    quadgk(λ-> cosh(β*γ*v*sqrt(ε^2-1)*sin(φ)*sin(X(ξ,ε,λ)))*sinh(β*γ*v*sqrt(ε^2-1)*cos(φ)*cos(X(ξ,ε,λ))),λ_c(ε),λ_max(ξ,ε))[1],ε_min(ξ), infinity)[1]
    return integral
end
function J_φ_ABS(ξ, φ)
   # coeff = -2 * α_sch * m_0^3*M / ξ 
    coeff = 1
    integral = coeff *quadgk(
        ε -> quadgk(
         λ-> λ/R(ξ, ε, λ) * exp(-β*γ*(ε+v*sqrt(ε^2-1)*cos(φ)*cos(X(ξ,ε,λ))))* 
         sinh(β*γ*v*sqrt(ε^2-1)*cos(φ)*cos(X(ξ,ε,λ))),0, λ_c(ε))[1],1,infinity)[1]
    return integral
end
function J_φ_SCATT(ξ, φ)
   # coeff = -4 * α_sch * m_0^3*M / ξ
    coeff = 1
    integral = coeff * quadgk(ε-> exp(-β*γ*ε) *  quadgk(λ-> λ* sinh(β*γ*v*sqrt(ε^2-1)*sin(φ)*sin(X(ξ,ε,λ))) * cosh(β*γ*v*sqrt(ε^2-1)*cos(φ)*cos(X(ξ,ε,λ))) / sqrt(R(ξ, ε, λ)),λ_c(ε),λ_max(ξ,ε))[1], ε_min(ξ), infinity)[1]
    return integral
end



