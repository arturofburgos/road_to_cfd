using Plots, LaTeXStrings, Measures # Possible to use PyPlot too.
include("schemes.jl")
include("plotdefaults.jl")
include("preprocessing.jl")
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
T = 3.0 # Total time
# Spatial domian number of points
nx = 201 # n + 1 = 201 -> n = 200

# Spatial grid for the numerical solution 
x = range(-1, 6, nx)
hx = x[2] - x[1]


# Define Mask:
# one way to ensure we never hit the boundary of the problem, 
# then in this specific case we have Dirichlet BCs since the 
# values are conserved throughout the whole time

# Create a boolean mask with all elements set to true
mask = trues(length(x))
# Set the first and last elements of the mask to false
# BCs
mask[1] = mask[end] = false


#-----------------------------------------------------------------------------#
# Note that we initialize the problem here in order to determine gamma
# We do that by defining the maximum fprime (in this case u) that
# we account when computing the CFL condition.
#-----------------------------------------------------------------------------#
condition = "shock"
u, u0 = initial_condition(x, condition)
fprime_max = maximum(abs.(fprime(u)))


# Define CFL
CFL = 0.95
ht = hx * CFL / (2 * fprime_max)
nt = round(T / ht)

# Print important info:
println("     T = $T")
println("tsteps = $nt")
println("    hx = $hx")
println("    ht = $ht")
println("   CFL = $CFL")
#-----------------------------------------------------------------------------------------------------------------------------------

#==================#
#    Processing:   #
#==================#
# Create Wave Shape:



# Plot initial condition
p = plot(x, u0, title="Initial Condition")
savefig(p, "$figures_dir/initial_condition.png")
display(p)


J = collect(1:nx)
Jm1 = circshift(J, 1)
Jp1 = circshift(J, -1)

println("J: ", J)
println("Jm1: ", Jm1)
println("Jp1: ", Jp1)

# Call update_solution (non and vectorized)
#u_hist = update_soln_vectorized2(nt, nx,  sigma)
u_hist, u_hist_exact = update_soln(nt, condition)


#-----------------------------------------------------------------------------------------------------------------------------------

#==================#
# Post Processing: #
#==================#
# Update plot function --> to be used in the .gif

# Using Multiple Dispatch:
function update_plot(n, u_hist, x)
    plot(x, u_hist[n])
end


function update_plot(n, u_hist, u_hist_exact, x)
    plot(x, u_hist[n])
    plot!(x, u_hist_exact[n])
end

# Create the .gif
anim1 = @animate for n in 1:length(u_hist)
    update_plot(n, u_hist, u_hist_exact, x)
end


# Save the .gif
gif(anim1, "$figures_dir/numerical_solution.gif", fps=15)
