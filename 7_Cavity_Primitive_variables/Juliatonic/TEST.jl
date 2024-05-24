include("FlowJulia.jl")
using .FlowJulia, Printf


length = 4.0
breadth = 4.0
colpts = 277
rowpts = 277
time = 50.0

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




CFL = CFL_number

SetTimeStep!(CFL, cavity)



setUBoundary!(cavity, noslip, noslip, flow, noslip)
setVBoundary!(cavity, noslip, noslip, noslip, noslip)
setPBoundary!(cavity, zeroflux, zeroflux, pressureatm, zeroflux)
GetStarredVelocities!(cavity, water)

SolvePressurePoisson!(cavity, water, zeroflux, zeroflux, pressureatm, zeroflux)
SolveMomentumEquation!(cavity, water)
SetCentrePUV!(cavity)



# print("end")
UpdateSolution!(; space=cavity, fluid=water, flow=flow, noslip=noslip, zeroflux=zeroflux, pressureatm=pressureatm,
    time=time, CFL_number=CFL_number, interval=interval, file_flag=file_flag)

