#=
Author: Ksymena Poradzisz
Updated: [2024-08-23]
Description: This is a code used to visualise a data obtain from  Schwarzschild.jl. 
=#
#=
To run code:
1. Update: in Linux shell:  juliaup update
2. Install missing packages:
	a) in Linux shell run julia
	b) in Julia REPL (shell) run: 
		import Pkg; Pkg.add("CSV", "PyCall", "DataFrames")
3. Run code from Linux shell: visualisation_schwarzschild.jl

Expected outcome:
A file named Schwarzschild_visualisation.mp4 saved in the path you run this code and the visualisation should pop in the window so after code is done running you would see the result. 
=#
using PyCall
using CSV, DataFrames

# Load the CSV file into a DataFrame
df = CSV.File("/home/korizekori/magisterka/Schwarzschild/data_Schwarzschild_2024-08-21_09-15-12.csv") |> DataFrame

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

# Visualisation in python
py"""
import numpy as np #numpy for mathematical functions such as sin, cos etc.
import matplotlib.pyplot as plt #for plotting
import matplotlib.animation as animation #for animations
from matplotlib.patches import Circle #for drawing black hole and orbits
import scipy.interpolate as sinter #for interpolating


def visualisation(x, y, J_x, J_y,n, beta, v): #definition of a function, which will create animation
    points = np.array(list(zip(x, y))) #this line create an array of points (x,y)
    ξ_hor = 2 # horizon
    ξ_ph = 3 #photon orbit
    
    # Create grids for streamplot - it is important to have points space evenly
    grid_x = np.linspace(min(x), max(x), 500)
    grid_y = np.linspace(min(y), max(y), 500)
    X, Y = np.meshgrid(grid_x, grid_y)

    # Perform interpolation to calculate values of Jx,Jy and n on meshgrid created above
    JX = sinter.griddata(points, J_x, (X, Y), method="nearest")
    JY = sinter.griddata(points, J_y, (X, Y), method="nearest")
    #particle density
    ni = sinter.griddata((x, y), n, (grid_x[None, :], grid_y[:, None]), method='cubic')
    # Function to return a value of interpolated vector field JX, JY on meshgrid X, Y for a given point (x0, y0)
    def vector_field(JX, JY, X, Y, x0, y0):
        ix = (np.abs(X[0] - x0)).argmin()  # Correct index for row
        iy = (np.abs(Y[:, 0] - y0)).argmin()  # Correct index for column
        return JX[iy, ix], JY[iy, ix] 
    # Function to create a streamline.
    #start_x and start_y are coordinates of the point in which the streamline will start
    #trail_length - specifies the maximum number of points in the particle's visible trail
    #(JX,JY) - vector field. Here - accretion current
    #X,Y - meshgrid
    #duration -parameter which specify duration of a line. Its unit is not exactly teremined as it is used to calculate number of frames and the real-time duration is determined by fps
    #dt - parameter which specify "time step" 
    def create_streamline(start_x, start_y, trail_length, duration,JX, JY, X, Y, dt):
        # Initial particle position
        px = np.array([start_x])
        py = np.array([start_y])

        # Lists to store trail positions
        trail_px = [px[0]]
        trail_py = [py[0]]

        frames = int(duration / dt)
        streamline_data = []

        for _ in range(frames):
            Vx, Vy = vector_field(JX, JY, X, Y, px, py)
            px += Vx * dt
            py += Vy * dt

            trail_px.append(px[0])
            trail_py.append(py[0])

            if len(trail_px) > trail_length:
                trail_px.pop(0)
                trail_py.pop(0)

            streamline_data.append((list(trail_px), list(trail_py)))

        return streamline_data
    #We want to have random initial positions of the streamlines, which are in the picture and which are in the range of calculated data. 
    def random_initial_position():
        x = np.random.uniform(-20, 20)
        y = np.random.uniform(-20, 20)
        return x,y
    # Set up the plot
    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    ax.set_xlim(-20, 20)
    ax.set_ylim(-20, 20)
    plt.title(r'$\beta = $' + f"{beta}"+ r', $v = $' + f"{v}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(True)

    #Showing particle density and creating legend for it (colorbar)
    cax = ax.imshow(ni, extent=[x.min(), x.max(), y.min(), y.max()], origin='lower', cmap='viridis', aspect='auto')
    colorbar = fig.colorbar(cax, ax=ax)
    colorbar.set_label(r'$\frac{n}{n_\infty}$', rotation=0, fontsize = 16,labelpad=15)

    # Initialize the plot elements: the lines representing the trails
    num_streamlines = 200
    lines = [ax.plot([], [], 'k-', lw=1, zorder = 1)[0] for _ in range(num_streamlines)]
    streamlines = [create_streamline(*random_initial_position(), trail_length=np.random.uniform(1000, 3000), duration=np.random.uniform(3, 10), JX = JX, JY = JY, X = X, Y = Y, dt =  0.01) for _ in range(num_streamlines)]
    streamline_frames = [0] * num_streamlines  # Tracks which frame each streamline is on

    # Function to initialize the animation
    def init():
        for line in lines:
            line.set_data([], [])
        return lines

    # Function to update the plot for each frame
    def update(frame):
        for i, (line, streamline_data) in enumerate(zip(lines, streamlines)):
            # Calculate the remaining lifespan of the streamline
            remaining_frames = len(streamline_data) - streamline_frames[i]

            # Fade out effect: gradually reduce the length or fade the color
            if remaining_frames < 120:  # Start fading out in the last 120 frames
                alpha = remaining_frames / 120.0
                line.set_alpha(alpha)
                
                # Gradually reduce the visible length of the trail
                if remaining_frames > 0:
                    trail_px, trail_py = streamline_data[streamline_frames[i]]
                    visible_length = max(1, int(len(trail_px) * remaining_frames / 120.0))
                    line.set_data(trail_px[-visible_length:], trail_py[-visible_length:])
                else:
                    line.set_data([], [])
            else:
                trail_px, trail_py = streamline_data[streamline_frames[i]]
                line.set_data(trail_px, trail_py)
                line.set_alpha(1.0)  # Reset alpha to fully visible

            streamline_frames[i] += 1

            # If the current streamline has expired, replace it with a new one
            if streamline_frames[i] >= len(streamline_data):
                streamlines[i] = create_streamline(*random_initial_position(), trail_length=np.random.uniform(1000, 3000), duration=np.random.uniform(3, 10), JX=JX, JY=JY, X=X, Y=Y, dt=0.01)
                streamline_frames[i] = 0

        return lines


    # Create the animation 
    ani = animation.FuncAnimation(fig, update, init_func=init, blit=True, interval=20, save_count=2000)
    # Add Black Hole with radius of ξ_hor and with center in (0,0)
    blackhole = Circle((0, 0), ξ_hor, edgecolor='black', facecolor='black',zorder = 100)
    ax.add_patch(blackhole)

    #Photon orbit
    Photon_orbit = Circle((0,0), ξ_ph, edgecolor = 'black', facecolor = 'none', linestyle = 'dotted')
    ax.add_patch(Photon_orbit) 

    #a vector showing in which direction the black hole is moving
    ax.quiver(0, 0, 5, 0, angles='xy', scale_units='xy', scale=1, color='black')

    # Save the animation
    ani.save('streamlines.mp4', writer='ffmpeg', fps=60, extra_args=['-vcodec', 'libx264', '-crf', '23'])

    # Display the animation (if needed in real-time)
    plt.show()

    
"""

py"visualisation"(x, y, J_X_TOTALsch, J_Y_TOTALsch,n, 2, -1/2)
