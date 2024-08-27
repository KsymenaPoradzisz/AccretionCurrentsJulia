using Symbolics, Dates
using DoubleExponentialFormulas, CSV
@variables M, m_0, ξ, ϵ_σ, ϵ_r, R̃, ε, U_λ, α,  β,v, λ, _ξ_, λ_c,φ

# Set global precision to 346 bits
setprecision(BigFloat, 426) # it corresponds to mathematica working precision 128
qde = QuadDE(BigFloat; h0=one(BigFloat)/8, maxlevel=10)
# Define the function X with high precision
X_and_Error(ξ, ε, λ, α) = qde(_ξ_ -> 
(BigFloat(λ) + (BigFloat(α) * BigFloat(ε)) / (BigFloat(1) - BigFloat(2) / BigFloat(_ξ_))) / 
(((BigFloat(α)^2) / (BigFloat(1) - BigFloat(2) / BigFloat(_ξ_)) + BigFloat(_ξ_)^2) * 
sqrt(BigFloat(ε)^2 - (BigFloat(1) - BigFloat(2) / BigFloat(_ξ_)) * 
(BigFloat(1) + BigFloat(λ)^2 / BigFloat(_ξ_)^2) - 
(BigFloat(α)^2 + BigFloat(2) * BigFloat(ε) * BigFloat(λ) * BigFloat(α)) / (BigFloat(_ξ_)^2))),
BigFloat(ξ), BigFloat(Inf); rtol = 1e-64
)



# Initialize data array
data = []

# Run calculations
for i in 1:1000
    # Generate random values with high precision
    ksi = abs(BigFloat(rand(9:100)))
    eps = abs(BigFloat(rand(2:20)))
    lambda = abs(BigFloat(rand(1:20)))
    alfa = abs(rand(BigFloat))
    #println(i, " and ", ksi)

    temp = []
    push!(temp, ksi)
    push!(temp, eps)
    push!(temp, lambda)
    push!(temp, alfa)  # Initial value for alpha

    try
        # Calculate the result and the elapsed time
        t = @elapsed X_and_Error(ksi, eps, lambda, 0)
        result,error = X_and_Error(ksi, eps, lambda, 0)


        if isinf(result)
            println("INF")
        else
            # Append results to the temp array
            push!(temp, result)
            push!(temp, t)
            #push!(temp, error)
            push!(data, Tuple(temp))
            i +=1
        end

        
    catch e_rror
        if isa(e_rror, DomainError)
            println("DomainError")
        end
    end
end

# Display the data
#println(data)
timestamp_for_file = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
CSV.write("/home/korizekori/magisterka/Julia_DATA_$(timestamp_for_file).csv", data, header = ["ξ","ε","λ", "α", "Result (Julia)", "Elapsed time in Julia[s]"])
