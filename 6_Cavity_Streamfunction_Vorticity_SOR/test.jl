include("function.jl")

using Plots, LaTeXStrings, Measures, .soln
#-----------------------------------------------------------------------------------------------------------------------------------
#=======================#
# Save files directory: #
#=======================#

# In order to save my plots define the main directory:
main_directory = dirname(@__FILE__)
# In case if we need to create the output file
mkpath("$main_directory/figures")
# Construct the path to the figures directory
figures_dir = joinpath(main_directory, "figures")
#-----------------------------------------------------------------------------------------------------------------------------------

#=======================#
# Initialize variables: #
#=======================#

# Define dimension of the cavity
L = 1

# Define the number of points of the domain
N = 31
nx = N
ny = N

# Spatial grid
x = range(0, L, N)
y = range(0, L, N)

h = x[2] - x[1]

# Properties of the simulation
U = 1
Re = 10
alpha = 1.5

# Initialize psi and w
psi = zeros(N, N)
w = zeros(N, N)

w_hist = []
psi_hist = []


# Update solution
begin
    psi_hist, w_hist, tn = update_soln(; nx, ny, psi, w, h, N, U, Re, alpha, w_hist, psi_hist)
end

#-----------------------------------------------------------------------------------------------------------------------------------
#=======================#
#    Post Processing:   #
#=======================#

# Compute velocities:
soln_size = length(psi_hist)

ux = []
uy = []


I = collect(1:nx)
Im1 = circshift(I, 1)
Ip1 = circshift(I, -1)

J = collect(1:ny)
Jm1 = circshift(J, 1)
Jp1 = circshift(J, -1)

for k = 1:soln_size
    ux_ij = (psi_hist[k][I, Jp1] - psi_hist[k][I, Jm1]) / (2 * h)
    uy_ij = -(psi_hist[k][Ip1, J] - psi_hist[k][Im1, J]) / (2 * h)


    push!(ux, ux_ij)
    push!(uy, uy_ij)

    # Update Velocity Boundary Conditions
    ux[k][:, end] .= U
    ux[k][:, 1] .= 0
    ux[k][1, :] .= 0
    ux[k][end, :] .= 0

    uy[k][:, end] .= 0
    uy[k][:, 1] .= 0
    uy[k][1, :] .= 0
    uy[k][end, :] .= 0
end



# Update plot functions --> to be used in the .gif
function update_vorticity_plot(n, w_hist, x, y)
    contourf(x, y, transpose(w_hist[n]),
        color=:jet, levels=20, aspect_ratio=:equal, clims=(-50.0, 15.0))
    xlabel!("x")
    ylabel!("y")
end


function update_streamfunction_plot(n, psi_hist, x, y)
    contourf(x, y, transpose(psi_hist[n]),
        color=:jet, levels=20, aspect_ratio=:equal, clims=(-0.12, 0.0))
    xlabel!("x")
    ylabel!("y")
end

function update_velocity_plot(n, ux, uy, x, y)

    # Create xp and yp as matrix as one the inputs of the quiver plot
    xp = x .+ 0 * y'
    yp = 0 * x .+ y'

    # Scale the ux and uy solution
    scale = 0.1
    ux_scaled = scale * ux[n]
    uy_scaled = scale * uy[n]
    quiver(xp, yp, quiver=(ux_scaled, uy_scaled), aspect_ratio=:equal, color=:black)
    xlabel!("x")
    ylabel!("y")
end

# Create the .gif's
anim1 = @animate for n in 1:10:length(w_hist)
    update_vorticity_plot(n, w_hist, x, y)
end

anim2 = @animate for n in 1:10:length(w_hist)
    update_streamfunction_plot(n, psi_hist, x, y)
end

anim3 = @animate for n in 1:10:length(w_hist)
    update_velocity_plot(n, ux, uy, x, y)
end


# Save the .gif
gif(anim1, "$figures_dir/vorticity.gif", fps=10)
gif(anim2, "$figures_dir/streamfunction.gif", fps=10)
gif(anim3, "$figures_dir/velocities.gif", fps=10)
