#=Author: Ksymena Poradzisz
Updated: 2024-08-30

Description: This is the file containing all necessary functions for calculating Kerr acrettion currents.
=#
using Symbolics
using QuadGK
using DoubleExponentialFormulas
using PolynomialRoots
using Polynomials
const infinity=123
#const infinity=Inf

println("Succesfully imported LIBKerr.jl")


using DoubleExponentialFormulas



# Define the function X with high precision

X_kerr(ξ, ε, λ, α) = quadde(_ξ_ -> 
(λ + (α * ε) / (1 - 2 / _ξ_)) / 
(((α^2) / (1 - 2 / _ξ_) + _ξ_^2) * 
sqrt(eps()+ε^2 - (1 - 2 / _ξ_) * (1 + λ^2 / _ξ_^2) - 
(α^2 + 2 * ε * λ * α) / _ξ_^2)),
ξ, Inf)[1]


precompile(X_kerr, (Float64, Float64, Float64, Float64))

U_λ_kerr(ξ, λ) = (1 - 2 / ξ) * (1 + λ^2 / ξ^2)
function ε_min_kerr(ξ, α, ϵ_σ)
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

function R̃_kerr(ξ, ε, α, λ, ϵ_σ)
    temp = ε^2 - U_λ_kerr(ξ, λ) - α * (α + 2 * ε * λ * ϵ_σ) / ξ^2
    if temp < 0
        return 0
    else
        return temp
    end
end

S(ξ, ε, λ, α, ϵ_σ, ϵ_r,φ) = exp(-(ε + ϵ_σ * v * sqrt(ε^2 - 1) * sin(φ - ϵ_σ * ϵ_r * (π / 2 - X_kerr(ξ, ε, λ, α)))) * β / sqrt(1 - v^2))

#=
function S(ξ, ε, λ, α, ϵ_σ, ϵ_r, φ)
    X = X_kerr(ξ, ε, λ, α)
   
    if isinf(X)
        X=0.0 # this is stupid solution...
    end
    
    return exp(-(ε + ϵ_σ * v * sqrt(ε^2 - 1) * sin(φ - ϵ_σ * ϵ_r * (π / 2 - X))) * β / sqrt(1 - v^2))
end
=#

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
    if R̃_kerr(ξ, ε, α, λ, ϵ_σ) <= 0
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
function J_t_ABS_kerr(ksi,φ, alfa, m_0)
    f = __jt_integrals__;
    temp1(λ, ε) = -m_0^3 / ksi * f(ksi, λ, ε, alfa,  1, -1,φ) #eps_sigma = 1; eps_r = -1
    temp2(λ, ε) = -m_0^3 / ksi * f(ksi, λ, ε, alfa, -1, -1,φ) #eps_sigma = -1; eps_r = -1
    result1, err1 = quadgk(ε -> quadgk(λ -> temp1(λ, ε), 0, λ_c_kerr(alfa, ε,  1, ksi))[1], 1, Inf)
    result2, err2 = quadgk(ε -> quadgk(λ -> temp2(λ, ε), 0, λ_c_kerr(alfa, ε, -1, ksi))[1], 1, Inf)
    result = result1 + result2
    return result

end
function J_φ_ABS_kerr(ksi,φ, alfa, m_0, M)
    f = __jφ_integrals__;
    temp1(λ, ε) = M * m_0^3 / ksi * f(ksi, λ, ε, alfa,  1, -1,φ) #eps_sigma = 1; eps_r = -1
    temp2(λ, ε) = M * m_0^3 / ksi * f(ksi, λ, ε, alfa, -1, -1,φ) #eps_sigma = -1; eps_r = -1
    result1, err1 = quadgk(ε -> quadgk(λ -> temp1(λ, ε), 0, λ_c_kerr(alfa, ε,  1, ksi))[1], 1, Inf)
    result2, err2 = quadgk(ε -> quadgk(λ -> temp2(λ, ε), 0, λ_c_kerr(alfa, ε, -1, ksi))[1], 1, Inf)
    result = result1 + result2
    return result
end
function J_r_ABS_kerr(ksi,φ, alfa, m_0)
    f = __jr_integrals__;
    temp1(λ, ε) = (m_0^3 * ksi) / (ksi * (ksi - 2) + alfa^2) * (f(ksi, λ, ε, alfa, 1, -1,φ)) #eps_sigma = 1; eps_r = -1
    temp2(λ, ε) = (m_0^3 * ksi) / (ksi * (ksi - 2) + alfa^2) * (f(ksi, λ, ε, alfa, -1, -1,φ)) #eps_sigma = -1; eps_r = -1
    result1, err = quadgk(ε -> quadgk(λ -> temp1(λ, ε), 0, λ_c_kerr(alfa, ε, 1, ksi))[1], 1, Inf)
    result2, err = quadgk(ε -> quadgk(λ -> temp2(λ, ε), 0, λ_c_kerr(alfa, ε, -1, ksi))[1], 1, Inf)
    result = result1 + result2
    return result
end

function J_t_SCATT_kerr(ksi,φ, alfa, m_0)
    f = __jt_integrals__;
    temp1(λ, ε) = -m_0^3 / ksi * (f(ksi, λ, ε, alfa, 1, 1,φ)) #eps_sigma = 1; eps_r = 1
    temp2(λ, ε) = -m_0^3 / ksi * (f(ksi, λ, ε, alfa, -1, 1,φ)) #eps_sigma = -1; eps_r = 1
    temp3(λ, ε) = -m_0^3 / ksi * (f(ksi, λ, ε, alfa, 1, -1,φ)) #eps_sigma = 1; eps_r = -1
    temp4(λ, ε) = -m_0^3 / ksi * (f(ksi, λ, ε, alfa, -1, -1,φ)) #eps_sigma = 1; eps_r = 1
    result1, err1 = quadgk(ε -> quadgk(λ -> temp1(λ, ε), λ_c_kerr(alfa, ε,  1, ksi), λ_max_kerr(ksi, ε, alfa,  1))[1], ε_min_kerr(ksi, alfa,  1), infinity) #lower boundary = λ_c; upper_boundary = λ_max 
    result2, err2 = quadgk(ε -> quadgk(λ -> temp2(λ, ε), λ_c_kerr(alfa, ε, -1, ksi), λ_max_kerr(ksi, ε, alfa, -1))[1], ε_min_kerr(ksi, alfa, -1), infinity)
    result3, err3 = quadgk(ε -> quadgk(λ -> temp3(λ, ε), λ_c_kerr(alfa, ε,  1, ksi), λ_max_kerr(ksi, ε, alfa,  1))[1], ε_min_kerr(ksi, alfa,  1), infinity)
    result4, err4 = quadgk(ε -> quadgk(λ -> temp4(λ, ε), λ_c_kerr(alfa, ε, -1, ksi), λ_max_kerr(ksi, ε, alfa, -1))[1], ε_min_kerr(ksi, alfa, -1), infinity)
    result = result1 + result2 + result3 + result4
    return result
end
function J_φ_SCATT_kerr(ksi,φ, alfa, m_0, M)
    f = __jφ_integrals__;
    temp1(λ, ε) = M * m_0^3 / ksi * (f(ksi, λ, ε, alfa, 1, 1,φ)) #eps_sigma = 1; eps_r = 1
    temp2(λ, ε) = M * m_0^3 / ksi * (f(ksi, λ, ε, alfa, -1, 1,φ)) #eps_sigma = -1; eps_r = 1
    temp3(λ, ε) = M * m_0^3 / ksi * (f(ksi, λ, ε, alfa, 1, -1,φ)) #eps_sigma = 1; eps_r = -1
    temp4(λ, ε) = M * m_0^3 / ksi * (f(ksi, λ, ε, alfa, -1, -1,φ)) #eps_sigma = 1; eps_r = 1t
    result1, err1 = quadgk(ε -> quadgk(λ -> temp1(λ, ε), λ_c_kerr(alfa, ε,  1, ksi), λ_max_kerr(ksi, ε, alfa,  1))[1], ε_min_kerr(ksi, alfa,  1),  infinity) #lower boundary = λ_c; upper_boundary = λ_max 
    result2, err2 = quadgk(ε -> quadgk(λ -> temp2(λ, ε), λ_c_kerr(alfa, ε, -1, ksi), λ_max_kerr(ksi, ε, alfa, -1))[1], ε_min_kerr(ksi, alfa, -1),  infinity)
    result3, err3 = quadgk(ε -> quadgk(λ -> temp3(λ, ε), λ_c_kerr(alfa, ε,  1, ksi), λ_max_kerr(ksi, ε, alfa,  1))[1], ε_min_kerr(ksi, alfa,  1),  infinity)
    result4, err4 = quadgk(ε -> quadgk(λ -> temp4(λ, ε), λ_c_kerr(alfa, ε, -1, ksi), λ_max_kerr(ksi, ε, alfa, -1))[1], ε_min_kerr(ksi, alfa, -1),  infinity)
    result = result1 + result2 + result3 + result4
    return result
end
function J_r_SCATT_kerr(ksi,φ, alfa, m_0)
    f = __jr_integrals__;
    temp1(λ, ε) = m_0^3 * ksi / (ksi * (ksi - 2) + alfa^2) * (f(ksi, λ, ε, alfa, 1, 1,φ)) #eps_sigma = 1; eps_r = 1
    temp2(λ, ε) = m_0^3 * ksi / (ksi * (ksi - 2) + alfa^2) * (f(ksi, λ, ε, alfa, -1, 1,φ)) #eps_sigma = -1; eps_r = 1
    temp3(λ, ε) = m_0^3 * ksi / (ksi * (ksi - 2) + alfa^2) * (f(ksi, λ, ε, alfa, 1, -1,φ)) #eps_sigma = 1; eps_r = -1
    temp4(λ, ε) = m_0^3 * ksi / (ksi * (ksi - 2) + alfa^2) * (f(ksi, λ, ε, alfa, -1, -1,φ)) #eps_sigma = 1; eps_r = 1
    result1, err1 = quadgk(ε -> quadgk(λ -> temp1(λ, ε), λ_c_kerr(alfa, ε,  1, ksi), λ_max_kerr(ksi, ε, alfa,  1))[1], ε_min_kerr(ksi, alfa,  1),  infinity) #lower boundary = λ_c; upper_boundary = λ_max 
    result2, err2 = quadgk(ε -> quadgk(λ -> temp2(λ, ε), λ_c_kerr(alfa, ε, -1, ksi), λ_max_kerr(ksi, ε, alfa, -1))[1], ε_min_kerr(ksi, alfa, -1),  infinity)
    result3, err3 = quadgk(ε -> quadgk(λ -> temp3(λ, ε), λ_c_kerr(alfa, ε,  1, ksi), λ_max_kerr(ksi, ε, alfa,  1))[1], ε_min_kerr(ksi, alfa,  1),  infinity)
    result4, err4 = quadgk(ε -> quadgk(λ -> temp4(λ, ε), λ_c_kerr(alfa, ε, -1, ksi), λ_max_kerr(ksi, ε, alfa, -1))[1], ε_min_kerr(ksi, alfa, -1),  infinity)
    result = result1 + result2 + result3 + result4
    return result
end
