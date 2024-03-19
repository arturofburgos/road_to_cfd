using Plots
#include("plotdefaults.jl")

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


#-----------------------------------------------------------------------------------------------------------------------------------

#=======================#
# Initialize variables: #
#=======================#
r_ini = 1.0; r_end = 2.0 # range in r-direction
θ_ini = 0.0; θ_end = pi/2 # range in θ-direction

a = 1 # Diffusivity
T = 10.0 # Total time
# Spatial domian number of points
nr = 31 # n + 1 = 31 -> n = 30
nθ = 31

# Spatial grid for the numerical solution 
r = range(r_ini, r_end, nr)
θ = range(θ_ini, θ_end, nθ)'
hr = r[2] - r[1]
hθ = θ[2] - θ[1]

# Temporal grid for the numerical solution


#-----------------------------------------------------------------------------------------------------------------------------------

#==================#
#    Processing:   #
#==================#
# Define BCs:

u_top = 90.0
u_bottom = 100.0
u_left = 90.0
u_right = 90.0

u = ones(nr, nθ)*90






# Create a filled contour plot
contourf(r, θ', transpose(u),
    title="Initial Condition", 
    color= :viridis)

display(gcf())





# Set labels for x and y axes


x1 = r.*cos.(θ) # =f(s,phi)
y1 = r.*sin.(θ) # =g(s,phi)


# Create a filled contour plot
contourf(x1, y1, transpose(u),
    title="Initial Condition", 
    color= :viridis)
display(gcf())

