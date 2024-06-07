timeinit = time()
include("FlowJulia.jl")
using .FlowJulia, Printf
using CairoMakie


length = 4.0
breadth = 4.0
colpts = 111
rowpts = 111
final_time = 50.0


#================================#
#       Define Miscelaneous      #
#================================#

CFL_number = 0.8 # Do not touch this unless solution diverges
file_flag = 1 # Keep 1 to print results to file
interval = 100 # Record values in file per interval number of iterations
plot_flag = 1 # Keep 1 to plot results at the end

#===============================#
#  Define Physical Parammeters  #
#===============================#

rho = 1.0
mu = 0.01

#================================#
#   Define Momentum Quantities   #
#================================#
u_in = 1.0
v_wall = 0.0
p_out = 0.0



#================================#
# Create Space and Fluid structs #
#================================#

cavity = Space(rowpts, colpts)
CreateMesh!(cavity, rowpts, colpts)
SetDeltas!(cavity, breadth, length)
water = Fluid(rho, mu)

#======================#
# Boundary Definitions #
#======================#

flow = Boundary("D", u_in)
noslip = Boundary("D", v_wall)
zeroflux = Boundary("N", 0.0)
pressureatm = Boundary("D", p_out)


# Begin the run:

println("######## Beginning FlowPy Simulation ########")
println("#############################################")
println("# Simulation time: $(final_time)")
println("# Mesh: $(colpts) x $(rowpts)")
println("# Re/u: $(rho * length / mu)\tRe/v: $(rho * breadth / mu)")
MakeResultDirectory()

# elapsed_time = @elapsed begin
#     UpdateSolution!(; space=cavity, fluid=water, flow=flow, noslip=noslip, zeroflux=zeroflux, pressureatm=pressureatm,
#     time=time, CFL_number=CFL_number, interval=interval, file_flag=file_flag)
# end

UpdateSolution!(; space=cavity, fluid=water, flow=flow, noslip=noslip, zeroflux=zeroflux, pressureatm=pressureatm, 
    time=final_time, CFL_number=CFL_number, interval=interval, file_flag=file_flag, timeinit=timeinit)


#=============================#
# Ghia results for comparison #
#=============================#

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


# Using Plots as a backend for plotting
# # Plot 1: Horizontal velocity vs Vertical distance along center
# p1 = plot(y, u_c[:, ceil(Int, colpts/2)], label="Numerical Solution", color=:darkblue, linewidth=2)
# scatter!(y_g, u_g, label="Ghia et al. Benchmark", color=:red, marker=:x, markersize=3)
# xlabel!("Vertical distance along center")
# ylabel!("Horizontal velocity")
# title!("Benchmark plot 1")

# # Plot 2: Vertical velocity vs Horizontal distance along center
# p2 = plot(x, v_c[ceil(Int, rowpts/2), :], label="Numerical Solution", color=:darkblue, linewidth=2)
# scatter!(x_g, v_g, label="Ghia et al. Benchmark", color=:red, marker=:x, markersize=3)
# xlabel!("Horizontal distance along center")
# ylabel!("Vertical velocity")
# title!("Benchmark plot 2")
# display(p1)
# display(p2)

# contourf(x,y, p_c, color=:viridis, levels=8, aspect_ratio=:equal, clims = (0.4, 2.0))
# xp = x .+ 0 * y'
# yp = 0 * x .+ y'
# scale = 0.1
# ux_scaled = scale * u_c
# vy_scaled = scale * v_c
# quiver!(xp, yp, quiver=(ux_scaled', vy_scaled'), aspect_ratio=:equal, color=:black)
# function meshgrid(vx::AbstractVector{T}, vy::AbstractVector{T}) where T
#     X = [x for _ in vy, x in vx]
#     Y = [y for y in vy, _ in vx]
#     X, Y
# end

# Using Makie as a backend for plotting
fig1 = Figure()
ax1 = Axis(fig1[1,1], 
          title = "Benchmark plot 1",
          xlabel = "Vertical distance along center",
          ylabel = "Horizontal velocity"
          )
lines!(ax1, y, u_c[:, ceil(Int, colpts/2)], label = "Numerical solution",color = :blue)
scatter!(ax1, y_g, u_g, marker = :utriangle, label = "Ghia et al. Benchmark", color = :red)
axislegend()
current_figure()
display(fig1)

fig2 = Figure()
ax2 = Axis(fig2[1,1], 
          title = "Benchmark plot 2",
          xlabel = "Horizontal distance along center",
          ylabel = "Vertical velocity"
          )
lines!(ax2, y, v_c[ceil(Int, rowpts/2), :], label = "Numerical solution",color = :blue)
scatter!(ax2, y_g, v_g, marker = :utriangle, label = "Ghia et al. Benchmark", color = :red)
axislegend()
current_figure()
display(fig2)


function VectorField(p ;length = length, breadth = breadth)
    # Scale p values to range from 0 to lenght_value and from 0 to breadth_value

    scaled_p1 = p[1] * colpts
    scaled_p2 = p[2] * rowpts
    # Clamp scaled values to ensure they stay within valid range
    i = clamp(round(Int, scaled_p1/length), 1, colpts)
    j = clamp(round(Int, scaled_p2/breadth), 1, rowpts)
    
    # Fetch the vector components from U and V matrices
    return Point2f(u_c[j, i], v_c[j, i])
end


fig3 = Figure()
ax3 = Axis(fig3[1,1], 
          title = "Pressure and Velocity Streamlines", 
          xlabel = "x",
          ylabel = "y"
          )
xlims!(ax3, 0.0, length)
ylims!(ax3, 0.0, breadth)

xinterval_vec_field = 0..length
yinterval_vec_field = 0..breadth

plt1 = contourf!(ax3, x, y, p_c', levels = 15)
plt2 = streamplot!(ax3, VectorField, xinterval_vec_field, yinterval_vec_field,
 colormap=:magma)
Colorbar(fig3[1, 2], plt1)
Colorbar(fig3[1, 3], plt2)
display(fig3)

fig4 = Figure()
ax4 = Axis(fig4[1,1], 
          title = "Pressure and Velocity Vectors", 
          xlabel = "x",
          ylabel = "y",
          )
xlims!(ax4, 0.0, length)
ylims!(ax4, 0.0, breadth)

xinterval = Int(round(rowpts/(length^2)))
yinterval = Int(round(colpts/(breadth^2)))
amplification_factor = 1

plt1 = contourf!(ax4, x, y, p_c', levels = 10)
plt2 = arrows!(ax4, x[1:xinterval:end], y[1:yinterval:end],
 (u_c[1:yinterval:end, 1:xinterval:end])'*amplification_factor,
 (v_c[1:yinterval:end, 1:xinterval:end])'*amplification_factor)
Colorbar(fig4[1, 2], plt1)
display(fig4)