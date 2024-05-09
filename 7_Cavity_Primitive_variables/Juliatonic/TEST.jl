include("FlowJulia.jl")
using .FlowJulia, Printf


length = 4.0
breadth = 4.0
colpts = 257
rowpts = 257
time = 1

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
p_out = 0.0

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

while t < time
    #printf("\rSimulation time left: %.2f", time - t)
    flush(stdout)

    CFL = CFL_number
    SetTimeStep!(CFL, cavity)
    timestep = cavity.dt

    setUBoundary!(cavity, noslip, noslip, flow, noslip)
    setVBoundary!(cavity, noslip, noslip, noslip, noslip)
    setPBoundary!(cavity, zeroflux, zeroflux, pressureatm, zeroflux)
    GetStarredVelocities!(cavity, water)

    SolvePressurePoisson!(cavity, water, zeroflux, zeroflux, pressureatm, zeroflux)
    SolveMomentumEquation!(cavity, water)

    SetCentrePUV!(cavity)
    # if file_flag == 1
    #     WriteToFile(cavity, i, interval)
    # end

    t += timestep
    # println(t)
    i += 1
    println(i)

    if i == 5
        println(cavity.u_c)
    elseif i == 55
        println(cavity.u_c)
    elseif i == 105
        println(cavity.u_c)
    elseif i == 155
        println(cavity.u_c)
    end
end
