#=
Author: Ksymena Poradzisz
Updated: [2024-09-10]
Description: This is a code used to visualise a data obtain from  Schwarzschild.jl. 
=#

println("Running visualisation with static streamlines")
using Plots, Glob
using FileIO
using CSV, PyCall
using DataFrames, Interpolations


directory = pwd()
pattern = "data_Schw_beta_*_v_*_dim_*_*.csv"
filelist = glob(pattern, directory)
if length(filelist) == 0
    error("No files found matching the pattern.")
end
# Sort the files by last modified time in descending order to get the newest file first
sorted_files = sort(filelist, by=file -> stat(file).mtime, rev=true)

# Select the newest file
filename = sorted_files[1]
println("Data file: $(filename)")


#read values beta and v from the filename
regex = r"beta_(\d+(\.\d+)?)_v_(-?\d+(\.\d+)?)_dim_(\d+(\.\d+)?)_"
match_result = match(regex, filename)
if match_result !== nothing
    beta_value = parse(Float64, match_result.captures[1])
    v_value = parse(Float64, match_result.captures[3])
    dim_value = parse(Float64, match_result.captures[5])
    println("Beta value extracted from filename: ", beta_value)
    println("V value extracted from filename: ", v_value)
    println("dim value extracted from filename: ", dim_value)
else
    error("Beta or V or dim value not found in the filename.")
end

df = CSV.File(filename) |> DataFrame

# Convert the columns to vectors
x = convert(Vector{Float64}, df.x)
y = convert(Vector{Float64}, df.y)
J_X_TOTAL_Schw = convert(Vector{Float64}, df.J_X_TOTAL_Schw)
J_Y_TOTAL_Schw = convert(Vector{Float64}, df.J_Y_TOTAL_Schw)
n = convert(Vector{Float64}, df.n)

#=
# Check for missing values and remove rows with missing values if any
if any(ismissing, x) || any(ismissing, y) || any(ismissing, J_X_TOTALsch) || any(ismissing, J_Y_TOTALsch) || any(ismissing, n)
    df = dropmissing(df)
    x = convert(Vector{Float64}, df.x)
    y = convert(Vector{Float64}, df.y)
    J_X_TOTALsch = convert(Vector{Float64}, df.J_X_TOTAlsch)
    J_Y_TOTALsch = convert(Vector{Float64}, df.J_Y_TOTALsch)
end
=#



# Visualisation in python
py"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as sinter
import matplotlib.animation as animation
from matplotlib.patches import Circle

def visualisation(x, y, J_x, J_y, n,beta,v,dim,save_path):
    # Create pair of points (x, y)
    points = list(zip(x, y))
    ξ_hor = 2 # horizon
    ξ_ph = 3 #photon orbit
    # Create grids for streamplot
    grid_x = np.linspace(min(x), max(x), 500)
    grid_y = np.linspace(min(y), max(y), 500)
    

    X, Y = np.meshgrid(grid_x, grid_y)
    
    ni = sinter.griddata((x, y), n, (grid_x[None, :], grid_y[:, None]), method='cubic')

    # Perform interpolation
    JX = sinter.griddata(points, J_x, (X, Y), method="cubic")
    JY = sinter.griddata(points, J_y, (X, Y), method="cubic")
    

    # Initialize the plot
    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    ax.set_xlim(-dim, dim)
    ax.set_ylim(-dim, dim)
    strm = ax.streamplot(X, Y, JX, JY, color='black', density=1, arrowsize=0, linewidth = 0.75,broken_streamlines=False)
    # Add Black Hole with radius of ξ_hor and with center in (0,0)
    blackhole = Circle((0, 0), ξ_hor, edgecolor='black', facecolor='black')
    ax.add_patch(blackhole)

    #Photon orbit
    Photon_orbit = Circle((0,0), ξ_ph, edgecolor = 'black', facecolor = 'none', linestyle = 'dotted')
    ax.add_patch(Photon_orbit)

    #a vector showing in which direction the black hole is moving
    ax.quiver(0, 0, 5, 0, angles='xy', scale_units='xy', scale=1, color='black')

    cax = ax.imshow(ni, extent=[x.min(), x.max(), y.min(), y.max()], origin='lower', cmap='viridis', aspect='auto')
    colorbar = fig.colorbar(cax, ax=ax)
    colorbar.set_label(r'$\frac{n}{n_\infty}$', rotation=0, fontsize = 16,labelpad=15)
    plt.title(r'$\beta = $' + f"{beta}"+ r', $v = $' + f"{v}")
    # Generate streamline paths
    paths = strm.lines.get_paths()
    segments = [path.vertices for path in paths]

    # Number of arrows
    num_arrows = 70
    
    # Initialize arrow positions and directions
    arrow_positions = np.zeros((num_arrows, 2))
    arrow_directions = np.zeros((num_arrows, 2))
    arrow_segments = [segments[np.random.choice(len(segments))] for _ in range(num_arrows)]
    arrow_offsets = np.random.rand(num_arrows)

    arrows = ax.quiver(arrow_positions[:, 0], arrow_positions[:, 1], 
                       arrow_directions[:, 0], arrow_directions[:, 1], 
                       angles='xy', scale_units='xy', scale=1.5, color='black')

    def update_arrows(num):
        for i in range(num_arrows):
            segment = arrow_segments[i]
            offset = arrow_offsets[i] + num * 0.001  # Move arrows along the path
            offset %= 1  # Wrap around

            # Get position along the path
            pos_idx = int(offset * (len(segment) - 1))
            pos = segment[pos_idx]
            if pos_idx + 1 < len(segment):
                next_pos = segment[pos_idx + 1]
            else:
                next_pos = segment[0]  # Wrap around to the start of the segment

            direction = next_pos - pos
            arrow_positions[i] = pos
            arrow_directions[i] = direction

        arrows.set_offsets(arrow_positions)
        arrows.set_UVC(arrow_directions[:, 0], arrow_directions[:, 1])
        return arrows,

    ani = animation.FuncAnimation(fig, update_arrows, frames=500, interval=10, blit=True)
    # Save animation as MP4
    ani.save(f"{save_path}.mp4", writer='ffmpeg', fps=60)
    
    # Save animation as GIF
    ani.save(f"{save_path}.gif", writer='pillow', fps=60)

   # plt.show() #uncomment if you want to see animation in real-time
"""

py"visualisation"(x, y, J_X_TOTAL_Schw, J_Y_TOTAL_Schw, n, beta_value, v_value, dim_value, "Schw_visualisation_static_$(beta_value)_$(v_value)_$(dim_value)")
