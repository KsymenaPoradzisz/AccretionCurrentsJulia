#=Author: Ksymena Poradzisz
Updated: 2024-09-09

Description: This is the file containing all necessary functions for calculating Schwarzschild acrettion currents.
=#

using Symbolics
using QuadGK
using DoubleExponentialFormulas
using PolynomialRoots
using Polynomials
const infinity = 123

println("Succesfully imported LIBSchw.jl")
#######SCHWARZSCHILD#############
λ_max_Schw(ξ, ε) = ξ * sqrt(ε^2 / (1 - 2 / ξ) - 1)
λ_c_Schw(ε) = sqrt(12 / (1 - 4 / (1 + 3 * ε / sqrt(9 * ε^2 - 8))^2))
U_λ_Schw(ξ, λ) = (1 - 2 / ξ) * (1 + λ^2 / ξ^2)
function R_Schw(ξ, ε, λ)
    temp = ε^2 - U_λ_Schw(ξ, λ)
    if temp <= 0
        return Inf
    else
        return temp
    end

end
X_Schw(ξ, ε, λ) = quadde(_ξ_-> λ/(_ξ_^2*sqrt(R_Schw(_ξ_, ε, λ))),ξ,Inf)[1]
function ε_min_Schw(ξ)
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

function J_t_ABS_Schw(ξ,φ)
    coeff = -2*M*m_0^3/ξ
  #  coeff = 1
    integral = coeff *  quadgk( ε-> ε* quadgk(λ-> cosh(β*γ*v*sqrt(ε^2-1)*sin(φ)*sin(X_Schw(ξ, ε, λ)))*
    exp(-β*γ*(ε+v*sqrt(ε^2-1)*cos(φ)*cos(X_Schw(ξ,ε,λ))))/sqrt(R_Schw(ξ, ε, λ)),0,λ_c_Schw(ε))[1],1,infinity)[1]
    return integral
end
function J_t_SCATT_Schw(ξ, φ)
    coeff = -4*M*m_0^3/ξ
    #coeff = 1
    integral=  coeff*quadgk( 
                ε-> exp(-β*γ*ε) * ε *
                 quadgk( λ ->  1*cosh(β*γ*v*sqrt(ε^2-1)*cos(φ)*cos(X_Schw(ξ, ε, λ))) *cosh(β*γ*v*sqrt(ε^2-1)*sin(φ)*sin(X_Schw(ξ, ε, λ)))/sqrt(R_Schw(ξ, ε, λ)), λ_c_Schw(ε), λ_max_Schw(ξ,ε))[1],ε_min_Schw(ξ), infinity)[1]
    return integral
end
function J_r_ABS_Schw(ξ,φ)
    coeff = -2*M*m_0^3/(ξ-2) 
    #coeff = 1
    integral = coeff *quadgk( ε->  
    quadgk(λ-> cosh(β*γ*v*sqrt(ε^2-1)*sin(φ)*sin(X_Schw(ξ, ε, λ)))*exp(-β*γ*(ε+v*sqrt(ε^2-1)*cos(φ)*cos(X_Schw(ξ,ε,λ)))),0,λ_c_Schw(ε))[1],1,infinity)[1]
    return integral
end
function J_r_SCATT_Schw(ξ,φ)
    coeff = 4*M*m_0^3/(ξ-2)
    #coeff = 1
    integral = coeff * quadgk(ε-> exp(-β*γ*ε) * 
    quadgk(λ-> cosh(β*γ*v*sqrt(ε^2-1)*sin(φ)*sin(X_Schw(ξ,ε,λ)))*sinh(β*γ*v*sqrt(ε^2-1)*cos(φ)*cos(X_Schw(ξ,ε,λ))),λ_c_Schw(ε),λ_max_Schw(ξ,ε))[1],ε_min_Schw(ξ), infinity)[1]
    return integral
end

function J_φ_ABS_Schw(ξ, φ)
    coeff = -2 * M * m_0^3*M / ξ 
   # coeff = 1
    integral = coeff *quadgk(
        ε -> quadgk(
         λ-> λ/sqrt(R_Schw(ξ, ε, λ)) * exp(-β*γ*(ε+v*sqrt(ε^2-1)*cos(φ)*cos(X_Schw(ξ,ε,λ))))* 
         sinh(β*γ*v*sqrt(ε^2-1)*sin(φ)*sin(X_Schw(ξ,ε,λ))),0, λ_c_Schw(ε))[1],1,infinity)[1]
    return integral
end
function J_φ_SCATT_Schw(ξ, φ)
    coeff = -4 * M * m_0^3*M / ξ
   # coeff = 1
    integral = coeff * quadgk(ε-> exp(-β*γ*ε) *  quadgk(λ-> λ* sinh(β*γ*v*sqrt(ε^2-1)*sin(φ)*sin(X_Schw(ξ,ε,λ))) * cosh(β*γ*v*sqrt(ε^2-1)*cos(φ)*cos(X_Schw(ξ,ε,λ))) / sqrt(R_Schw(ξ, ε, λ)),λ_c_Schw(ε),λ_max_Schw(ξ,ε))[1], ε_min_Schw(ξ), infinity)[1]
    return integral
end
