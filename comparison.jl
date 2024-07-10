#=
Author: Ksymena Poradzisz
Updated: [2024-07-10]
=#


using Symbolics
using QuadGK
using DoubleExponentialFormulas
using PolynomialRoots
using Polynomials

const v = -0.5 # predkosc w jedn. c
const β = 2  # Thermodynamic β=1/(kT), aka coolness
const γ = 1 / sqrt(1 - v^2) #global usage unnecessary
const infinity = 123
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
    coeff = -2*M*m_0^3/ξ
  #  coeff = 1
    integral = coeff *  quadgk( ε-> ε* quadgk(λ-> cosh(β*γ*v*sqrt(ε^2-1)*sin(φ)*sin(X_sch(ξ, ε, λ)))*
    exp(-β*γ*(ε+v*sqrt(ε^2-1)*cos(φ)*cos(X_sch(ξ,ε,λ))))/sqrt(R_sch(ξ, ε, λ)),0,λ_c_sch(ε))[1],1,infinity)[1]
    return integral
end
function J_t_SCATT_sch(ξ, φ)
    coeff = -4*M*m_0^3/ξ
    #coeff = 1
    integral=  coeff*quadgk( 
                ε-> exp(-β*γ*ε) * ε *
                 quadgk( λ ->  1*cosh(β*γ*v*sqrt(ε^2-1)*cos(φ)*cos(X_sch(ξ, ε, λ))) *cosh(β*γ*v*sqrt(ε^2-1)*sin(φ)*sin(X_sch(ξ, ε, λ)))/sqrt(R_sch(ξ, ε, λ)), λ_c_sch(ε), λ_max_sch(ξ,ε))[1],ε_min_sch(ξ), infinity)[1]
    return integral
end
function J_r_ABS_sch(ξ,φ)
    coeff = -2*M*m_0^3/(ξ-2) 
    #coeff = 1
    integral = coeff *quadgk( ε->  
    quadgk(λ-> cosh(β*γ*v*sqrt(ε^2-1)*sin(φ)*sin(X_sch(ξ, ε, λ)))*exp(-β*γ*(ε+v*sqrt(ε^2-1)*cos(φ)*cos(X_sch(ξ,ε,λ)))),0,λ_c_sch(ε))[1],1,infinity)[1]
    return integral
end
function J_r_SCATT_sch(ξ,φ)
    coeff = -4*M*m_0^3/(ξ-2)
    #coeff = 1
    integral = coeff * quadgk(ε-> exp(-β*γ*ε) * 
    quadgk(λ-> cosh(β*γ*v*sqrt(ε^2-1)*sin(φ)*sin(X_sch(ξ,ε,λ)))*sinh(β*γ*v*sqrt(ε^2-1)*cos(φ)*cos(X_sch(ξ,ε,λ))),λ_c_sch(ε),λ_max_sch(ξ,ε))[1],ε_min_sch(ξ), infinity)[1]
    return integral
end

function J_φ_ABS_sch(ξ, φ)
    coeff = -2 * M * m_0^3*M / ξ 
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

#############KERR################
X_kerr(ξ, ε, λ, α) = quadde(_ξ_ -> (λ + (α * ε) / (1 - 2 / _ξ_)) / (((α^2) / (1 - 2 / _ξ_) + _ξ_^2) * sqrt((ε^2 - (1 - 2 / _ξ_) * (1 + λ^2 / _ξ_^2) - (α^2 + 2 * ε * λ * α) / (_ξ_^2)))), ξ, Inf)[1]

U_λ_kerr(ξ, λ) = (1 - 2 / ξ) * (1 + λ^2 / ξ^2)


function R̃_kerr(ξ, ε, α, λ, ϵ_σ)
    temp = ε^2 - U_λ_kerr(ξ, λ) - α * (α + 2 * ε * λ * ϵ_σ) / ξ^2
    if temp < 0
        return 0
    else
        return temp
    end
end

S(ξ, ε, λ, α, ϵ_σ, ϵ_r,φ) = exp(-(ε + ϵ_σ * v * sqrt(ε^2 - 1) * sin(φ - ϵ_σ * ϵ_r * (π / 2 - X_kerr(ξ, ε, λ, α)))) * β / sqrt(1 - v^2))

function λ_max_kerr(ξ, ε, α, ϵ_σ)
    if ξ == 2 && sign(α) * ϵ_σ == 1
        temp = -(α^2 - 4 * ε^2) / (2 * ε * α)
    else
        temp = ξ / (ξ - 2) * (-ϵ_σ * ε * α + sqrt(α^2 * (ε^2 + 2 / ξ - 1) + (ξ - 2) * (ξ * (ε^2 - 1) + 2)))
    end
    return temp
end

function λ_c_kerr(α, ε, ϵ_σ, ξ)
    # println("alfa = ", α, " ε =  ", ε, "ϵ_σ = " ,ϵ_σ, " ξ = ", ξ)
    poly_coeffs = [-α^4 - α^6 * (-1 + ε^2),
        (-4 * ϵ_σ * α^3 * ε - 6 * ϵ_σ * α^5 * ε * (-1 + ε^2)),
        (16 - 2 * α^2 - 4 * ϵ_σ^2 * α^2 * ε^2 + 18 * α^2 * (-1 + ε^2) - 3 * α^4 * (-1 + ε^2) - 12 * ϵ_σ^2 * α^4 * ε^2 * (-1 + ε^2)),
        (-4 * ϵ_σ * α * ε + 36 * ϵ_σ * α * ε * (-1 + ε^2) - 12 * ϵ_σ * α^3 * ε * (-1 + ε^2) - 8 * ϵ_σ^3 * α^3 * ε^3 * (-1 + ε^2)),
        (-1 + 18 * (-1 + ε^2) - 3 * α^2 * (-1 + ε^2) - 12 * ϵ_σ^2 * α^2 * ε^2 * (-1 + ε^2) + 27 * (-1 + ε^2)^2),
        6 * ϵ_σ * α * ε * (-1 + ε^2),
        1 - ε^2]
    poly = Polynomial(poly_coeffs)
    sols_imaginary = PolynomialRoots.roots(Float64.(poly_coeffs))
    sols = filter(sol -> abs(imag(sol)) < 1e-15, sols_imaginary)
    limit_λ = -ϵ_σ * α + 2 + 2 * sqrt(1 - ϵ_σ * α)
    if isempty(sols)
        return nothing
    else
        real_sols = broadcast(abs, (real.(sols)))
        maks = maximum(real_sols)
        return maks >= limit_λ ? maks : nothing
    end
end

function __jt_integrals__(ξ, λ, ε, α, ϵ_σ, ϵ_r,φ)
    if R̃_kerr(ξ, ε, α, λ, ϵ_σ) == 0
        return 0
    else
        temp = ε * S(ξ, ε, λ, α, ϵ_σ, ϵ_r,φ) / sqrt(R̃_kerr(ξ, ε, α, λ, ϵ_σ))
        return temp
    end
end
function __jr_integrals__(ξ, λ, ε, α, ϵ_σ, ϵ_r,φ)
    if R̃_kerr(ξ, ε, α, λ, ϵ_σ) == 0
        return 0
    else
        temp = ϵ_r * S(ξ, ε, λ, α, ϵ_σ, ϵ_r,φ)
        return temp
    end
end
function __jφ_integrals__(ξ, λ, ε, α, ϵ_σ, ϵ_r,φ)
    R = R̃_kerr(ξ, ε, α, λ, ϵ_σ)
    
    if R̃_kerr(ξ, ε, α, λ, ϵ_σ) == 0
        return 0
    else
        temp = ϵ_σ * S(ξ, ε, λ, α, ϵ_σ, ϵ_r,φ) * (λ + ϵ_σ * α * ε) / sqrt(R̃_kerr(ξ, ε, α, λ, ϵ_σ))
        #println("R = $(R) for ξ = $(ξ), λ = $(λ), ε = $(ε), ϵ_σ = $(ϵ_σ), funkcja podcałkowa = $(temp)")
        return temp
    end
end
function J_t_ABS_kerr(f, ksi,φ, alfa, m_0)
    temp1(λ, ε) = -m_0^3 / ksi * f(ksi, λ, ε, alfa,  1, -1,φ) #eps_sigma = 1; eps_r = -1
    temp2(λ, ε) = -m_0^3 / ksi * f(ksi, λ, ε, alfa, -1, -1,φ) #eps_sigma = -1; eps_r = -1
    result1, err1 = quadgk(ε -> quadgk(λ -> temp1(λ, ε), 0, λ_c_kerr(alfa, ε,  1, ksi))[1], 1, Inf)
    result2, err2 = quadgk(ε -> quadgk(λ -> temp2(λ, ε), 0, λ_c_kerr(alfa, ε, -1, ksi))[1], 1, Inf)
    result = result1 + result2
    return result

end
function J_φ_ABS_kerr(f, ksi,φ, alfa, m_0, M)
    temp1(λ, ε) = M * m_0^3 / ksi * f(ksi, λ, ε, alfa,  1, -1,φ) #eps_sigma = 1; eps_r = -1
    temp2(λ, ε) = M * m_0^3 / ksi * f(ksi, λ, ε, alfa, -1, -1,φ) #eps_sigma = -1; eps_r = -1
    result1, err1 = quadgk(ε -> quadgk(λ -> temp1(λ, ε), 0, λ_c_kerr(alfa, ε,  1, ksi))[1], 1, Inf)
    result2, err2 = quadgk(ε -> quadgk(λ -> temp2(λ, ε), 0, λ_c_kerr(alfa, ε, -1, ksi))[1], 1, Inf)
    result = result1 + result2
    return result
end
function J_r_ABS_kerr(f, ksi,φ, alfa, m_0)
    temp1(λ, ε) = (m_0^3 * ksi) / (ksi * (ksi - 2) + alfa^2) * (f(ksi, λ, ε, alfa, 1, -1,φ)) #eps_sigma = 1; eps_r = -1
    temp2(λ, ε) = (m_0^3 * ksi) / (ksi * (ksi - 2) + alfa^2) * (f(ksi, λ, ε, alfa, -1, -1,φ)) #eps_sigma = -1; eps_r = -1
    result1, err = quadgk(ε -> quadgk(λ -> temp1(λ, ε), 0, λ_c_kerr(alfa, ε, 1, ksi))[1], 1, Inf)
    result2, err = quadgk(ε -> quadgk(λ -> temp2(λ, ε), 0, λ_c_kerr(alfa, ε, -1, ksi))[1], 1, Inf)
    result = result1 + result2
    return result
end


for x in 12:-0.5:5 # printing differences between two integrals 

    J_t_ABS_Diff = J_t_ABS_kerr(__jt_integrals__,x,1,0,1)   - J_t_ABS_sch(x,1)
    J_r_ABS_Diff = J_r_ABS_kerr(__jr_integrals__,x,1,0,1)   - J_r_ABS_sch(x,1)
    J_φ_ABS_Diff = J_φ_ABS_kerr(__jφ_integrals__,x,1,0,1,1) - J_φ_ABS_sch(x,1)

    J_t_ABS_Rel = J_t_ABS_kerr(__jt_integrals__,x,1,0,1)   / J_t_ABS_sch(x,1)
    J_r_ABS_Rel = J_r_ABS_kerr(__jr_integrals__,x,1,0,1)   / J_r_ABS_sch(x,1)
    J_φ_ABS_Rel = J_φ_ABS_kerr(__jφ_integrals__,x,1,0,1,1) / J_φ_ABS_sch(x,1)

    println("x = $x , J_t_ABS: Kerr / Sch = ", J_t_ABS_Rel)
    println("x = $x , J_r_ABS: Kerr / Sch = ", J_r_ABS_Rel)
    println("x = $x , J_φ_ABS: Kerr / Sch = ", J_φ_ABS_Rel)    
    println("")

end
