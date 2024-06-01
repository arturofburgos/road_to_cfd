include("FlowJulia.jl")
using .FlowJulia, Printf#
using CairoMakie


length = 1.0
breadth = 1.0
colpts = 31
rowpts = 31
time = 500.0

###############################MISC############################################
CFL_number = 0.8 # Do not touch this unless solution diverges
file_flag = 1 # Keep 1 to print results to file
interval = 100 # Record values in file per interval number of iterations
plot_flag = 1 # Keep 1 to plot results at the end

###########################DEFINE PHYSICAL PARAMETERS##########################
rho = 1.0
mu = 0.01

##########################DEFINE INITIAL MOMENTUM PARAMETERS###################
u_in = 1.0
v_wall = 0.0
p_out = 1.0

###############################################################################
########################CREATE SPACE OBJECT####################################
cavity = Space(rowpts, colpts)
CreateMesh!(cavity, rowpts, colpts)
SetDeltas!(cavity, breadth, length)
water = Fluid(rho, mu)

###############################################################################
#########################BOUNDARY DEFINITIONS##################################
########################CREATE BOUNDARY OBJECTS################################
###########################VELOCITY############################################
flow = Boundary("D", u_in)
noslip = Boundary("D", v_wall)
zeroflux = Boundary("N", 0.0)
############################PRESSURE###########################################
pressureatm = Boundary("D", p_out)

#######################USER INPUT ENDS#########################################
###############################################################################
#############################INITIALIZATION####################################
t = 0
i = 0
############################THE RUN############################################
println("######## Beginning FlowPy Simulation ########")
println("#############################################")
println("# Simulation time: $(time)")
println("# Mesh: $(colpts) x $(rowpts)")
println("# Re/u: $(rho * length / mu)\tRe/v: $(rho * breadth / mu)")
#println("# Save outputs to text file: $(bool(file_flag))")
MakeResultDirectory()


function SetTimeStep!(CFL::Float64, space::Space)
    dt = CFL / sum([maximum(space.u) / space.dx, maximum(space.v) / space.dy])
    # Escape condition if dt is infinity due to zero velocity initially
    if isinf(dt)
        dt = CFL * (space.dx + space.dy)
    end
    space.dt = dt
end

# print("end")
UpdateSolution!(; space=cavity, fluid=water, flow=flow, noslip=noslip, zeroflux=zeroflux, pressureatm=pressureatm,
    time=time, CFL_number=CFL_number, interval=interval, file_flag=file_flag)




x = range(0, length , colpts)
y = range(0, breadth, rowpts)

u = cavity.u
v = cavity.v
p = cavity.p
u_c = cavity.u_c
v_c = cavity.v_c
p_c = cavity.p_c


y_g = [0, 0.0547, 0.0625, 0.0703, 0.1016, 0.1719, 0.2813, 0.4531, 0.5, 0.6172, 0.7344, 0.8516, 0.9531, 0.9609, 0.9688, 0.9766]
u_g = [0, -0.08186, -0.09266, -0.10338, -0.14612, -0.24299, -0.32726, -0.17119, -0.11477, 0.02135, 0.16256, 0.29093, 0.55892, 0.61756, 0.68439, 0.75837]

x_g = [0, 0.0625, 0.0703, 0.0781, 0.0983, 0.1563, 0.2266, 0.2344, 0.5, 0.8047, 0.8594, 0.9063, 0.9453, 0.9531, 0.9609, 0.9688]
v_g = [0, 0.1836, 0.19713, 0.20920, 0.22965, 0.28124, 0.30203, 0.30174, 0.05186, -0.38598, -0.44993, -0.23827, -0.22847, -0.19254, -0.15663, -0.12146]

# Scale the benchmark coordinates
y_g = [breadth * y for y in y_g]
x_g = [length * x for x in x_g]



# Plot 1: Horizontal velocity vs Vertical distance along center
p1 = plot(y, u_c[:, ceil(Int, colpts/2)], label="Numerical Solution", color=:darkblue, linewidth=2)
scatter!(y_g, u_g, label="Ghia et al. Benchmark", color=:red, marker=:x, markersize=3)
xlabel!("Vertical distance along center")
ylabel!("Horizontal velocity")
title!("Benchmark plot 1")


# Plot 2: Vertical velocity vs Horizontal distance along center
p2 = plot(x, v_c[ceil(Int, rowpts/2), :], label="Numerical Solution", color=:darkblue, linewidth=2)
scatter!(x_g, v_g, label="Ghia et al. Benchmark", color=:red, marker=:x, markersize=3)
xlabel!("Horizontal distance along center")
ylabel!("Vertical velocity")
title!("Benchmark plot 2")

display(p1)
display(p2)



# contourf(x,y, p_c, color=:viridis, levels=8, aspect_ratio=:equal, clims = (0.4, 2.0))

# xp = x .+ 0 * y'
# yp = 0 * x .+ y'
# scale = 0.1

# ux_scaled = scale * u_c

# vy_scaled = scale * v_c

# quiver!(xp, yp, quiver=(ux_scaled', vy_scaled'), aspect_ratio=:equal, color=:black)



function meshgrid(vx::AbstractVector{T}, vy::AbstractVector{T}) where T
    X = [x for _ in vy, x in vx]
    Y = [y for y in vy, _ in vx]
    X, Y
end



function vector_field5(p)
    # Convert the coordinates to indices
    i = clamp(round(Int, p[1]*31), 1, 31)  # X index
    j = clamp(round(Int, p[2]*31), 1, 31)  # Y index
    
    # Fetch the vector components from U and V matrices
    return Point2f(u_c[j, i], v_c[j, i])
end


# vector_field5([31,0])


# fig, ax, plt = streamplot(vector_field5, -0..31, -0..31, colormap=:magma)

# fig, ax, plt = streamplot(vector_field5, -0..1, -0..1, colormap=:magma)


# streamplot!(vector_field5, -0..31, -0..31, colormap=:magma)






fig, ax, plt = contourf(x, y, p_c', levels=25)
streamplot!(vector_field5, 0..1, 0..1, colormap=:magma)
Colorbar(fig[1, 2], plt)

fig