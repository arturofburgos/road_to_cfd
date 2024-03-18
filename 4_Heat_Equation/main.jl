using Plots, LaTeXStrings, Measures # Possible to use PyPlot too.
include("schemes.jl")
include("plotdefaults.jl")

#When we uncomment the lline below when generating the plots we
#ended up getting this warning: 
#┌Warning: Multiple series with different fill alpha share a colorbar. 
#Colorbar may not reflect all series correctly.
#In this case we leave it commented.
# gr() # one of the Plots BackEnd


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
Lx = 1;
Ly = 1;
a = 1 # Diffusivity
T = 10.0 # Total time
# Spatial domian number of points
nx = 31 # n + 1 = 31 -> n = 30
ny = 31

# Spatial grid for the numerical solution 
x = range(0, Lx, nx)
y = range(0, Ly, ny)
hx = x[2] - x[1]
hy = y[2] - y[1]

# Temporal grid for the numerical solution
ht = (1 / (2 * a)) * (((hx * hy)^2) / (hx^2 + hy^2))
nt = round(T / ht)

# Print important info:
println("     T = $T")
println("tsteps = $nt")
println("    hx = $hx")
println("    hy = $hy")
println("    ht = $ht")
#-----------------------------------------------------------------------------------------------------------------------------------

#==================#
#    Processing:   #
#==================#
# Define BCs:

u_top = 50.0
u_bottom = 50.0
u_left = 0.0
u_right = 100.0

u0 = ones(nx, ny)*0
u = apply_BC(u0)


# Plot initial condition:
fig = plot(size=(550, 400))

# Set labels for x and y axes
xlabel!("x")
ylabel!("y")

# Create a filled contour plot
contourf!(x, y, transpose(u),
    title="Initial Condition", 
    color= :viridis)
savefig(fig, "$figures_dir/initial_condition.png")
display(fig)


I = collect(1:nx)
Im1 = circshift(I, 1)
Ip1 = circshift(I, -1)

J = collect(1:ny)
Jm1 = circshift(J, 1)
Jp1 = circshift(J, -1)


println("  I: ", I)
println("Im1: ", Im1)
println("Ip1: ", Ip1)
println("  J: ", J)
println("Jm1: ", Jm1)
println("Jp1: ", Jp1)

# Call update_solution
u_hist= update_soln(nt)


#-----------------------------------------------------------------------------------------------------------------------------------

#==================#
# Post Processing: #
#==================#
# Update plot function --> to be used in the .gif
function update_plot(n, u_hist, x, y)
    contourf(x, y, transpose(u_hist[n]), 
    color= :viridis)
    xlabel!("x")
    ylabel!("y")
end


# Create the .gif
anim1 = @animate for n in 1:500
    update_plot(n, u_hist, x, y)
end


# Save the .gif
gif(anim1, "$figures_dir/numerical_solution.gif", fps=30)
