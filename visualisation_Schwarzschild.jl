#=
Author: Ksymena Poradzisz
Updated: [2024-08-09]
Description: This is a code used to visualise a data obtain from  Schwarzschild.jl. 
=#
#=
To run code:
1. Update: in Linux shell:  juliaup update
2. Install missing packages:
	a) in Linux shell run julia
	b) in Julia REPL (shell) run: 
		import Pkg; Pkg.add("Plots", "CSV", "PyCall", "DataFrames", "Interpolations", "Glob")
3. Run code from Linux shell: visualisation_schwarzschild.jl

Expected outcome:
A file named Schwarzschild_visualisation.mp4 saved in the path you run this code and the visualisation should pop in the window so after code is done running you would see the result. 
=#
using Plots, Glob
using CSV, PyCall
using DataFrames, Interpolations

@pyimport matplotlib.pyplot as plt
@pyimport numpy as np
@pyimport scipy.interpolate as sinter
@pyimport matplotlib.animation as animation

# Load the CSV file into a DataFrame
directory = pwd()
pattern = "data_Schwarzschild_beta_*_v_*_*.csv"
filelist = glob(pattern, directory)
if length(filelist) == 0
    error("No files found matching the pattern.")
end
filename = filelist[1]
println("Data file: $(filename)")

df = CSV.File(filename) |> DataFrame

# Convert the columns to vectors
x = convert(Vector{Float64}, df.x)
y = convert(Vector{Float64}, df.y)
J_X_TOTALsch = convert(Vector{Float64}, df.J_X_TOTAlsch)
J_Y_TOTALsch = convert(Vector{Float64}, df.J_Y_TOTALsch)
n = convert(Vector{Float64}, df.n)

# Check for missing values and remove rows with missing values if any
if any(ismissing, x) || any(ismissing, y) || any(ismissing, J_X_TOTALsch) || any(ismissing, J_Y_TOTALsch) || any(ismissing, n)
    df = dropmissing(df)
    x = convert(Vector{Float64}, df.x)
    y = convert(Vector{Float64}, df.y)
    J_X_TOTALsch = convert(Vector{Float64}, df.J_X_TOTAlsch)
    J_Y_TOTALsch = convert(Vector{Float64}, df.J_Y_TOTALsch)
end

#read values beta and v from the filename
regex = r"beta_(\d+)_v_(-?\d+\.\d+)_"
match_result = match(regex, filename)
if match_result !== nothing
    beta_value = match_result.captures[1]
    v_value = match_result.captures[2]
    println("Beta value extracted from filename: ", beta_value)
    println("V value extracted from filename: ", v_value)
else
    error("Beta or V value not found in the filename.")
end


# Visualisation in python
py"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as sinter
import matplotlib.animation as animation
from matplotlib.patches import Circle

def visualisation(x, y, J_x, J_y, n,beta,v,save_path):
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
    ax.set_xlim(-20, 20)
    ax.set_ylim(-20, 20)
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
    
  #  try:
   #     if save_format == 'mp4':
   #         ani.save(save_path, writer='ffmpeg', fps=60)
   #     elif save_format == 'gif':
   #         ani.save(save_path, writer='pillow', fps=60)
   # except ValueError as e:
   #     print(f"Error saving animation: {e}. Trying to save as GIF instead.")
  #      ani.save(save_path.replace('.mp4', '.gif'), writer='pillow', fps=30)

    plt.show()
"""

py"visualisation"(x, y, J_X_TOTALsch, J_Y_TOTALsch, n, beta_value, v_value, "Schwarzschild_visualisation_$(beta_value)_$(v_value)")
