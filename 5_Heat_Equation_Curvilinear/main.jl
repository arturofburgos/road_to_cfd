using Plots, LaTeXStrings, Measures # Possible to use PyPlot too.
pyplot() # PyPlot  BackEnd
include("schemes.jl")

#-----------------------------------------------------------------------------------------------------------------------------------
#=======================#
# Save files directory: #
#=======================#

# In order to save my plots define the main directory:
main_directory = dirname(@__FILE__)
# Construct the path to the figures directory
figures_dir = joinpath(main_directory, "figures")
#-----------------------------------------------------------------------------------------------------------------------------------

#=======================#
# Initialize variables: #
#=======================#
r_ini = 1.0; r_end = 2.0; # range in r-direction
θ_ini = 0.0; θ_end = pi / 2; # range in θ-direction

a = 1 # Diffusivity
T = 10.0 # Total time
# Spatial domian number of points
nr_true = 31
nθ_true = 31

nr = nr_true+2 # +2 because of the ghost cells
nθ = nθ_true+2 # +2 because of the ghost cells


# Spatial grid for the numerical solution 
r = range(r_ini, r_end, nr) 
θ = range(θ_ini, θ_end, nθ)' # we  transpose it in order to create a meshgrid when initializing x and y
hr = r[2] - r[1]
hθ = θ[2] - θ[1]



# Define Mask:
# one way to ensure we never hit the boundary of the problem, 
# then in this specific case we have Dirichlet BCs since the 
# values are conserved throughout the whole time

# Create a boolean mask with all elements set to false
mask = falses(length(r),length(θ))
# Set the inner elements of the mask to true
# Inner elements
mask[2:nr-1, 2:nθ-1] .= true



# Temporal grid for the numerical solution
ht = (1 / (2 * a)) * (((hr * hθ)^2) / (hr^2 + hθ^2))
nt = round(T / ht)

# Print important info:
println("     T = $T")
println("tsteps = $nt")
println("    hr = $hr")
println("    hθ = $hθ")
println("    ht = $ht")
#-----------------------------------------------------------------------------------------------------------------------------------

#==================#
#    Processing:   #
#==================#
# Define BCs:

u_top = 90.0
u_bottom = 100.0
u_left = 90.0
u_right = 90.0

u0 = ones(nr, nθ) * 90
u = apply_BC(u0)




# Plot initial condition in POLAR DOMAIN Coordinates:
fig1 = plot(size=(550, 400))

# Set labels for x and y axes
xlabel!("r")
ylabel!("θ")

# Create a filled contour plot
contourf!(r, θ, transpose(u),
    title="Initial Condition",
    color=:jet)
savefig(fig1, "$figures_dir/initial_condition_polar.png")
display(fig1)




# Plot initial condition in PHYSICAL DOMAIN Coordinates:
fig2 = plot(size=(550, 400))

# Set labels for x and y axes
xlabel!("x")
ylabel!("y")

# Physical Domain Coordinates
x = r .* cos.(θ) # NxN
y = r .* sin.(θ) # NxN

# Important to notice that, in this case, contourf with the BackEnd of pyplot
# allows to get both types of inputs, if it is polar, we get a vector for r and θ
# if it is physical coordinates it gets a matrix NxN for x and y, similar to
# np.meshgrid in python

# Create a filled contour plot
contourf!(x, y, transpose(u),
    title="Initial Condition",
    color=:jet, levels=20, aspect_ratio=:equal, xlim = (0,2), ylim = (0,2))
savefig(fig2, "$figures_dir/initial_condition_physical.png")
display(fig2)




# Defining vecotrized indices
I = collect(1:nr)
Im1 = circshift(I, 1)
Ip1 = circshift(I, -1)

J = collect(1:nθ)
Jm1 = circshift(J, 1)
Jp1 = circshift(J, -1)


println("  I: ", I)
println("Im1: ", Im1)
println("Ip1: ", Ip1)
println("  J: ", J)
println("Jm1: ", Jm1)
println("Jp1: ", Jp1)

#-----------------------------------------------------------------------------------------------------------------------------------

#==================#
# Update Solution: #
#==================#

# Call update_solution
u_hist = update_soln(nt)

#-----------------------------------------------------------------------------------------------------------------------------------

#==================#
# Post Processing: #
#==================#
# Update plot functions --> to be used in the .gif
function update_plot_physical(n, u_hist, x, y)
    contourf(x, y, transpose(u_hist[n]),
        color=:jet, levels=20, aspect_ratio=:equal)
    xlabel!("x")
    ylabel!("y")
end

function update_plot_polar(n, u_hist, x, y)
    contourf(x, y, transpose(u_hist[n]),
        color=:jet, levels=20)
    xlabel!("r")
    ylabel!("θ")
end


# Create the .gif's
anim1 = @animate for n in 1:100:length(u_hist)
    update_plot_physical(n, u_hist, x, y)
end

anim2 = @animate for n in 1:100:length(u_hist)
    update_plot_polar(n, u_hist, r, θ)
end


# Save the .gif
gif(anim1, "$figures_dir/numerical_solution_physical_domain.gif", fps=10)
gif(anim2, "$figures_dir/numerical_solution_polar_domain.gif", fps=10)

