using Symbolics
using Dates, CSV
using DataFrames
using CSV
using CairoMakie
using DataFrames
using LaTeXStrings
using StatsBase  , Glob, FileIO


#search for files of data to compare
directory = pwd()
pattern_Schw = "DATA_Schw_comp_*.csv"
pattern_kerr = "DATA_Kerr_comp_*.csv"

Schwlist = glob(pattern_Schw, directory)
if length(Schwlist) == 0
    error("No files found matching the Schwarzschild pattern.")
end
Kerrlist = glob(pattern_kerr, directory)
if length(Schwlist) == 0
    error("No files found matching the Kerr pattern.")
end

sorted_Schw = sort(Schwlist, by=file -> stat(file).mtime, rev=true)
sorted_Kerr = sort(Kerrlist, by=file -> stat(file).mtime, rev=true)

file_Schw = sorted_Schw[1]
file_Kerr = sorted_Kerr[1]


data_Schw = CSV.File(file_Schw) |> DataFrame
data_Kerr = CSV.File(file_Kerr) |> DataFrame

#checking whether the r and phi columns are the same, because otherwise comparison makes no sense

tol = 1e-9 

r_is_equal = all(abs.(data_Schw.r .- data_Kerr.r) .< tol)
φ_is_equal = all(abs.(data_Schw.φ  .- data_Kerr.φ ) .< tol)


if r_is_equal  &&  φ_is_equal 
    println("The columns r and φ are identical. Comparing integrals...")
else
    error("The columns r and φ are not identical. Can't proceed")
end

#uploading data from files
for name in names(data_Schw)
   # println(name)
    @eval $(Symbol(name)) = convert(Vector{Float64}, data_Schw.$name)
end

for name in names(data_Kerr)
    #println(name)
    @eval $(Symbol(name)) = convert(Vector{Float64}, data_Kerr.$name)
end


# Compute the values of relative J_*

J_t_ABS_Rel = ifelse.(abs.(J_t_ABSkerr) .< tol .|| abs.(J_t_ABS_Schw) .< tol, NaN, J_t_ABSkerr ./ J_t_ABS_Schw)
J_r_ABS_Rel = ifelse.(abs.(J_r_ABSkerr) .< tol .|| abs.(J_r_ABS_Schw) .< tol, NaN, J_r_ABSkerr ./ J_r_ABS_Schw)
J_φ_ABS_Rel = ifelse.(abs.(J_φ_ABSkerr) .< tol .|| abs.(J_φ_ABS_Schw) .< tol, NaN, J_φ_ABSkerr ./ J_φ_ABS_Schw)

J_t_SCATT_Rel = ifelse.(abs.(J_t_SCATTkerr) .< tol .|| abs.(J_t_SCATT_Schw) .< tol, NaN, J_t_SCATTkerr ./ J_t_SCATT_Schw)
J_r_SCATT_Rel = ifelse.(abs.(J_r_SCATTkerr) .< tol .|| abs.(J_r_SCATT_Schw) .< tol, NaN, J_r_SCATTkerr ./ J_r_SCATT_Schw)
J_φ_SCATT_Rel = ifelse.(abs.(J_φ_SCATTkerr) .< tol .|| abs.(J_φ_SCATT_Schw) .< tol, NaN, J_φ_SCATTkerr ./ J_φ_SCATT_Schw)


#calculating difference and its log10 
ABS_t_list = ifelse.(isnan.(J_t_ABS_Rel), NaN, log10.(abs.(J_t_ABS_Rel .- 1)))
ABS_r_list = ifelse.(isnan.(J_r_ABS_Rel), NaN, log10.(abs.(J_r_ABS_Rel .- 1)))
ABS_φ_list = ifelse.(isnan.(J_φ_ABS_Rel), NaN, log10.(abs.(J_φ_ABS_Rel .- 1)))

SCATT_t_list = ifelse.(isnan.(J_t_SCATT_Rel), NaN, log10.(abs.(J_t_SCATT_Rel .- 1)))
SCATT_r_list = ifelse.(isnan.(J_r_SCATT_Rel), NaN, log10.(abs.(J_r_SCATT_Rel .- 1)))
SCATT_φ_list = ifelse.(isnan.(J_φ_SCATT_Rel), NaN, log10.(abs.(J_φ_SCATT_Rel .- 1)))



df_transformed = DataFrame(
    r = r, 
    φ = φ, 
    ABS_t = ABS_t_list,
    ABS_r = ABS_r_list,
    ABS_φ = ABS_φ_list,
    SCATT_t = SCATT_t_list,
    SCATT_r = SCATT_r_list,
    SCATT_φ = SCATT_φ_list
)



#Transforming data to make it easy to draw
pivot_t_abs = unstack(df_transformed, :r, :φ, :ABS_t)
pivot_r_abs = unstack(df_transformed, :r, :φ, :ABS_r)
pivot_φ_abs = unstack(df_transformed, :r, :φ, :ABS_φ)
pivot_t_scatt = unstack(df_transformed, :r, :φ, :SCATT_t)
pivot_r_scatt = unstack(df_transformed, :r, :φ, :SCATT_r)
pivot_φ_scatt = unstack(df_transformed, :r, :φ, :SCATT_φ)

ABS_t_matrix = Matrix(pivot_t_abs[:, 2:end])
ABS_r_matrix = Matrix(pivot_r_abs[:, 2:end])
ABS_φ_matrix = Matrix(pivot_φ_abs[:, 2:end])
SCATT_t_matrix = Matrix(pivot_t_scatt[:, 2:end])
SCATT_r_matrix = Matrix(pivot_r_scatt[:, 2:end])
SCATT_φ_matrix = Matrix(pivot_φ_scatt[:, 2:end])
#function allowing to save 2D histograms 
function save_heatmap(matrix, title, filename)
    xticks = sort(unique(r))
    yticks = sort(unique(φ))
    fig = Figure()
    yticks_labels = [round(y/π, digits=2) for y in yticks]
    ax = Axis(fig[1, 1], title = title, xlabel = "r", ylabel = "φ",xticks = (1:length(xticks), string.(xticks)), yticks = (1:length(yticks), string.(yticks_labels).*"π"))
    hm = heatmap!(ax, 1:length(xticks), 1:length(yticks), matrix, colormap = :viridis)
    Colorbar(fig[1, 2], hm, label = "Log(Difference from 1)")
    save(filename, fig)
end
#creating histograms and saving them
save_heatmap(ABS_t_matrix, L"J^{(abs)}_{t,relative}", "ABS_t.png")
save_heatmap(ABS_r_matrix, L"J^{(abs)}_{r,relative}", "ABS_r.png")
save_heatmap(ABS_φ_matrix, L"J^{(abs)}_{\phi,relative}", "ABS_φ.png")

save_heatmap(SCATT_t_matrix, L"J^{(scatt)}_{t,relative}", "SCATT_t.png")
save_heatmap(SCATT_r_matrix, L"J^{(scatt)}_{r,relative}", "SCATT_r.png")
save_heatmap(SCATT_φ_matrix, L"J^{(scatt)}_{\phi,relative}", "SCATT_φ.png")
