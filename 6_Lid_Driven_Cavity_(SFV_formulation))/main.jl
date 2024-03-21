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
a = 0.1 # Diffusivity
T = 10.0 # Total time
# Spatial domian number of points
nx = 41 # n + 1 = 31 -> n = 30
ny = 41
n = 41

# Spatial grid for the numerical solution 
x = range(0, Lx, nx)
y = range(0, Ly, ny)
# hx = x[2] - x[1]
# hy = y[2] - y[1]
h = x[2] - x[1]



# Temporal grid for the numerical solution
# ht = (1 / (2 * a)) * (((hx * hy)^2) / (hx^2 + hy^2))
ht = (1 / (2 * a)) * (((h * h)^2) / (h^2 + h^2))
nt = round(T / ht)

# Print important info:
println("     T = $T")
println("tsteps = $nt")
println("     h = $h")
# println("    hx = $hx")
# println("    hy = $hy")
println("    ht = $ht")


U = 1
Re = 10
alpha = 1.5
tn = 0 
psi = zeros(length(x),length(y))
w = zeros(length(x),length(y))
#-----------------------------------------------------------------------------------------------------------------------------------
# k = ones(length(x),length(y)).*2
# o = ones(length(x),length(y)).*1.99999999999
# L1norm(k,o)




#==================#
#    Processing:   #
#==================#
# Define BCs:

psi_top = (-2/h^2)*psi[2:end-1,end-1] .- 2*U/h
psi_bottom = (-2/h^2)*psi[2:end-1,2]
psi_left = (-2/h^2)*psi[2,2:end-1]
psi_right = (-2/h^2)*psi[end-1,2:end-1]

psi0 = ones(nx, ny)*0
w = apply_BC(psi)




# Plot initial condition for vorticity (w):
fig = plot(size=(550, 400))

# Set labels for x and y axes
xlabel!("x")
ylabel!("y")

# Create a filled contour plot
contourf!(x, y, transpose(w),
    title="Initial Condition", 
    color= :viridis)
savefig(fig, "$figures_dir/initial_condition.png")
display(fig)




# Plot initial condition for streamfunction (psi):
fig = plot(size=(550, 400))

# Set labels for x and y axes
xlabel!("x")
ylabel!("y")

# Create a filled contour plot
contourf!(x, y, transpose(psi),
    title="Initial Condition", 
    color= :viridis)
savefig(fig, "$figures_dir/initial_condition.png")
display(fig)


l1_target = 0.9
l1_norm_w = 1
l1_norm_psi = 1


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
w, psi = update_soln(psi, w, tn)


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
