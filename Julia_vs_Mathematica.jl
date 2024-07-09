using Plots
using LinearAlgebra
using Plots.PlotMeasures
using CSV
using Dates
using DataFrames

# Load the CSV file into a DataFrame
df = CSV.File("/home/korizekori/magisterka/Julia_vs_Mathematica20240709.csv") |> DataFrame
current_date = Dates.format(Dates.now(), "yyyy-mm-dd")
# Access the columns
time_rel = 1.0 ./df[!, "Time: Julia / Mathematica"]
result_rel = df[!, "Result: Julia / Mathematica"]

# Create a histogram which shows relative time:  Time in Mathematica / Time in Julia
hist1 = histogram(time_rel, title="Time:  Mathematica/Julia ", xlabel="How many times faster Julia is", ylabel="Counts", bins=20,legend=false)
save_path1 = "/home/korizekori/time_rel_$(current_date).png"
#println("Histogram saved to: ", save_path)
# Save the histogram to a specific location
savefig(hist1, save_path1)


# Create a histogram which shows relative result: Result in Mathematica/Result in Julia
hist2 = histogram(result_rel, title="Result:  Julia/Mathematica ", xlabel="Ratio of results", ylabel="Counts", bins=21,legend=false)
save_path2 = "/home/korizekori/result_rel_$(current_date).png"
#println("Histogram saved to: ", save_path)
# Save the histogram to a specific location
savefig(hist2, save_path2)
