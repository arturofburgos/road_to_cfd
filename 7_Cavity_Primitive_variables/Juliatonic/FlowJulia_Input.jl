include("FlowJulia.jl")
using .FlowJulia


length = 4.0  # Length of computational domain in the x-direction
breadth = 4.0  # Breadth of computational domain in the y-direction
colpts = 10  # Number of grid points in the x-direction #KEEP ODD
rowpts = 10  # Number of grid points in the y-direction #KEEP ODD


rho = 1.0
mu = 0.01
space = Space(rowpts, colpts)
CreateMesh!(space, rowpts, colpts)
SetDeltas!(space, breadth, length)
fluid = Fluid(rho, mu)


SetTimeStep!(0.8, space)

GetStarredVelocities!(space, fluid)
u_in=1.0 #Lid velocity
v_wall=0.0 #Velocity of fluid at the walls
p_out=1.0
flow=Boundary("D",u_in)
noslip=Boundary("D",v_wall)
zeroflux=Boundary("N",0.0)
pressureatm=Boundary("D",p_out)

setPBoundary!(space, zeroflux, zeroflux, pressureatm, zeroflux)


SolvePressurePoisson!(space, fluid, zeroflux, zeroflux, pressureatm, zeroflux)