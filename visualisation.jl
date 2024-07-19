#Visualisation of Kerr data
using Plots
using LinearAlgebra
using Plots.PlotMeasures
using CSV
using DataFrames
df = CSV.File("/home/korizekori/magisterka/Kerr/result_beta=2_v=-0.5_alpha=0.0001.csv") |> DataFrame
x = convert(Vector{Float64}, df.x)
y = convert(Vector{Float64}, df.y)
J_X_ABS = convert(Vector{Float64}, df.J_X_ABS)
J_Y_ABS = convert(Vector{Float64}, df.J_Y_ABS)
# Check for missing values
if any(ismissing, x) || any(ismissing, y) || any(ismissing, J_X_ABS) || any(ismissing, J_Y_ABS)
    # Remove rows with missing values
    df = dropmissing(df)
    x = convert(Vector{Float64}, df.x)
    y = convert(Vector{Float64}, df.y)
    J_X_ABS = convert(Vector{Float64}, df.J_X_ABS)
    J_Y_ABS = convert(Vector{Float64}, df.J_Y_ABS)
end

# Create the quiver plot with consistent arrow styles
Plots.quiver(x, y, quiver=(J_X_ABS, J_Y_ABS), quiver_scale=0.5, arrow=:arrow, linecolor=:red)


# Define circle parameters
circle_radius = 5
theta = range(0, stop=2π, length=100)
circle_x = circle_radius * cos.(theta)
circle_y = circle_radius * sin.(theta)

# Add the black full circle
Plots.plot!(circle_x, circle_y, linecolor=:black, linewidth=2, fillcolor=:black, fillalpha=1, label="", ratio=:equal)
title!("Vector Field Visualization")
Plots.xlabel!("x")
Plots.ylabel!("y")
