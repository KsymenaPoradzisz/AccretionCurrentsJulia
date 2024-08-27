# Define a type alias to toggle precision
const MyFloat = BigFloat  # Use BigFloat for high precision
#const MyFloat = Float64  # Use Float64 for machine precision


using Dates
using DoubleExponentialFormulas, CSV


# Check if MyFloat is BigFloat and set precision accordingly
if MyFloat === BigFloat
    setprecision(MyFloat, 426)  # Set precision to 426 bits (128 decimal digits) for BigFloat
end

qde = QuadDE(MyFloat; h0=one(MyFloat)/8, maxlevel=10)
# Define the function X with high precision
X_and_Error(ξ, ε, λ, α) = qde(_ξ_ -> 
(MyFloat(λ) + (MyFloat(α) * MyFloat(ε)) / (MyFloat(1) - MyFloat(2) / MyFloat(_ξ_))) / 
(((MyFloat(α)^2) / (MyFloat(1) - MyFloat(2) / MyFloat(_ξ_)) + MyFloat(_ξ_)^2) * 
sqrt(MyFloat(ε)^2 - (MyFloat(1) - MyFloat(2) / MyFloat(_ξ_)) * 
(MyFloat(1) + MyFloat(λ)^2 / MyFloat(_ξ_)^2) - 
(MyFloat(α)^2 + MyFloat(2) * MyFloat(ε) * MyFloat(λ) * MyFloat(α)) / (MyFloat(_ξ_)^2))),
MyFloat(ξ), MyFloat(Inf); rtol = 1e-64
)



# Initialize data array
data = []

# Run calculations
for i in 1:10000
    # Generate random values with high precision
    ksi    = MyFloat(rand(9:100))
    eps    = MyFloat(rand(2:20))
    lambda = MyFloat(rand(1:20))
    alfa   = MyFloat(rand(-15:15))/MyFloat(16)

    temp = []
    push!(temp, ksi)
    push!(temp, eps)
    push!(temp, alfa)  # Initial value for alpha

    try
        # Calculate the result and the elapsed time
        t = @elapsed ((result, error) = X_and_Error(ksi, eps, lambda, alfa))


        if isinf(result)
            println("INF")
        else
            # Append results to the temp array
            push!(temp, result)
            push!(temp, t)
            push!(data, Tuple(temp))
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
CSV.write("Julia_DATA_$(timestamp_for_file).csv", data, header = ["ξ","ε","λ", "α", "Result (Julia)", "Elapsed time in Julia[s]"], encoding="UTF-8")