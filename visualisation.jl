using Plots
using CSV, PyCall
using DataFrames, Interpolations

@pyimport matplotlib.pyplot as plt
@pyimport numpy as np
@pyimport scipy.interpolate as sinter
@pyimport matplotlib.animation as animation

# Load the CSV file into a DataFrame
df = CSV.File("/home/korizekori/magisterka/Schwarzschild/data_Schwarzschild_2024-07-29_10-33-13.csv") |> DataFrame

# Convert the columns to vectors
x = convert(Vector{Float64}, df.x)
y = convert(Vector{Float64}, df.y)
J_X_TOTALsch = convert(Vector{Float64}, df.J_X_TOTAlsch)
J_Y_TOTALsch = convert(Vector{Float64}, df.J_Y_TOTALsch)

# Check for missing values and remove rows with missing values if any
if any(ismissing, x) || any(ismissing, y) || any(ismissing, J_X_TOTALsch) || any(ismissing, J_Y_TOTALsch)
    df = dropmissing(df)
    x = convert(Vector{Float64}, df.x)
    y = convert(Vector{Float64}, df.y)
    J_X_TOTALsch = convert(Vector{Float64}, df.J_X_TOTAlsch)
    J_Y_TOTALsch = convert(Vector{Float64}, df.J_Y_TOTALsch)
end

# Define the Python visualization function
py"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as sinter
import matplotlib.animation as animation

def visualisation(x, y, J_x, J_y, save_path, save_format='mp4'):
    # Create pair of points (x, y)
    points = list(zip(x, y))
    
    # Create grids for streamplot
    grid_x = np.linspace(min(x), max(x), 280)
    grid_y = np.linspace(min(y), max(y), 280)
    
    # Create meshgrid
    X, Y = np.meshgrid(grid_x, grid_y)
    
    # Perform interpolation
    JX = sinter.griddata(points, J_x, (X, Y), method="cubic")
    JY = sinter.griddata(points, J_y, (X, Y), method="cubic")
    
    # Initialize the plot
    fig, ax = plt.subplots()
    strm = ax.streamplot(X, Y, JX, JY, color='black', density=1.5, arrowsize=0)

    # Generate streamline paths
    paths = strm.lines.get_paths()
    segments = [path.vertices for path in paths]
    
    # Number of arrows
    num_arrows = 100
    
    # Initialize arrow positions and directions
    arrow_positions = np.zeros((num_arrows, 2))
    arrow_directions = np.zeros((num_arrows, 2))
    arrow_segments = [segments[np.random.choice(len(segments))] for _ in range(num_arrows)]
    arrow_offsets = np.random.rand(num_arrows)

    # Create quiver for arrows
    arrows = ax.quiver(arrow_positions[:, 0], arrow_positions[:, 1], 
                       arrow_directions[:, 0], arrow_directions[:, 1], 
                       angles='xy', scale_units='xy', scale=1, color='black')

    def update_arrows(num):
        for i in range(num_arrows):
            segment = arrow_segments[i]
            offset = arrow_offsets[i] + num * 0.01  # Move arrows along the path
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

    ani = animation.FuncAnimation(fig, update_arrows, frames=200, interval=30, blit=True)

    try:
        if save_format == 'mp4':
            ani.save(save_path, writer='ffmpeg', fps=30)
        elif save_format == 'gif':
            ani.save(save_path, writer='imagemagick', fps=30)
    except ValueError as e:
        print(f"Error saving animation: {e}. Trying to save as GIF instead.")
        ani.save(save_path.replace('.mp4', '.gif'), writer='imagemagick', fps=30)

    plt.show()
"""

# Call the Python visualization function with data from Julia
py"visualisation"(x, y, J_X_TOTALsch, J_Y_TOTALsch, "output.mp4", "mp4")
