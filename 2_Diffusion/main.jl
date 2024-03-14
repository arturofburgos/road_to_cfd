using Plots, LaTeXStrings, Measures # Possible to use PyPlot too.
include("schemes.jl")
include("plotdefaults.jl")
gr()

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
a = 0.01 # Diffusivity
T = 1.0 # Total time
# Spatial domian number of points
nx = 65 # n + 1 = 65 -> n = 64

# Spatial grid for the numerical solution 
x = range(0, 1, nx)
hx = x[2] - x[1]


# Define sigma
sigma = 0.4
ht = (hx^2) * sigma / a
nt = round(T / ht)

# Print important info:
println("     a = $a")
println("     T = $T")
println("tsteps = $nt")
println("    hx = $hx")
println("    ht = $ht")
println(" sigma = $sigma")
#-----------------------------------------------------------------------------------------------------------------------------------

#==================#
#    Processing:   #
#==================#
# Create Wave Shape:
function square(x)
    u = zeros(size(x))
    u[findall(x -> x > 0.4 && x < 0.6, x)] .= 2
    return u
end

# Define the initial condition for numerical solution:
f = square
uinit = f(x)


# Plot initial condition
p = plot(x, uinit, title="Initial Condition")
savefig(p, joinpath(figures_dir, "initial_condition.png"))
display(p)

# In case if we want to use the vectorized functions 
# (and periodic boundary conditions), we define:
J = collect(1:nx)
Jm1 = circshift(J, 1)
Jp1 = circshift(J, -1)

println("J: ", J)
println("Jm1: ", Jm1)
println("Jp1: ", Jp1)

# Call update_solution (non and vectorized)
#u_hist = update_soln_vectorized2(nt, nx,  sigma)
u_hist = update_soln_vectorized1(nt, sigma)


#-----------------------------------------------------------------------------------------------------------------------------------

#==================#
# Post Processing: #
#==================#
# Update plot function --> to be used in the .gif

# Using Multiple Dispatch:
function update_plot(n, u_hist, x)
    plot(x, u_hist[n])
end


# Create the .gif
anim1 = @animate for n in 1:100
    update_plot(n, u_hist, x)
end


# Save the .gif
gif(anim1, joinpath(figures_dir, "dirichlet.gif"), fps=15)
