using CSV
using CairoMakie
using DataFrames
using LaTeXStrings
using StatsBase  # For unique sorting of ticks

# Load the DataFrame from the CSV file
df = CSV.File("/home/korizekori/magisterka/comparison_data_2024-08-05_10-45-36.csv") |> DataFrame

# Convert DataFrame columns to Vector{Float64}
r = convert(Vector{Float64}, df.r)
φ = convert(Vector{Float64}, df.φ)
J_t_ABS_rel = convert(Vector{Float64}, df.J_t_ABS_rel)
J_r_ABS_rel = convert(Vector{Float64}, df.J_r_ABS_Rel)
J_φ_ABS_rel = convert(Vector{Float64}, df.J_φ_ABS_Rel)

# Transform the values into log10(abs(value - 1))
ABS_t_list = log10.(abs.(J_t_ABS_rel .- 1))
ABS_r_list = log10.(abs.(J_r_ABS_rel .- 1))
ABS_φ_list = log10.(abs.(J_φ_ABS_rel .- 1))

# Create a DataFrame with the transformed values and coordinates
df_transformed = DataFrame(r = r, φ = φ, ABS_t = ABS_t_list, ABS_r = ABS_r_list, ABS_φ = ABS_φ_list)

# Define xticks and yticks based on unique sorted values of r and φ
xticks = sort(unique(r))
yticks = sort(unique(φ))

# Create pivot tables for each ABS list
pivot_t = unstack(df_transformed, :r, :φ, :ABS_t)
pivot_r = unstack(df_transformed, :r, :φ, :ABS_r)
pivot_φ = unstack(df_transformed, :r, :φ, :ABS_φ)

# Extract the matrices from the pivot tables
ABS_t_matrix = Matrix(pivot_t[:, 2:end])
ABS_r_matrix = Matrix(pivot_r[:, 2:end])
ABS_φ_matrix = Matrix(pivot_φ[:, 2:end])


function save_heatmap(matrix, title, filename)
    fig = Figure()
    yticks_labels = [round(y/π, digits=2) for y in yticks]
    ax = Axis(fig[1, 1], title = title, xlabel = "r", ylabel = "φ",xticks = (1:length(xticks), string.(xticks)), yticks = (1:length(yticks), string.(yticks_labels).*"π"))
    hm = heatmap!(ax, 1:length(xticks), 1:length(yticks), matrix, colormap = :viridis)
    Colorbar(fig[1, 2], hm, label = "Log(Difference from 1)")
    save(filename, fig)
end
save_heatmap(ABS_t_matrix, L"J^{(abs)}_{t,relative}", "ABS_t.png")
save_heatmap(ABS_r_matrix, L"J^{(abs)}_{r,relative}", "ABS_r.png")
save_heatmap(ABS_φ_matrix, L"J^{(abs)}_{φ,relative}", "ABS_φ.png")