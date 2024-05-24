include("FlowJulia.jl")
using .FlowJulia


length = 4.0  # Length of computational domain in the x-direction
breadth = 4.0  # Breadth of computational domain in the y-direction
colpts = 10  # Number of grid points in the x-direction #KEEP ODD
rowpts = 10  # Number of grid points in the y-direction #KEEP ODD
CFL_number = 0.8

rho = 1.0
mu = 0.01
space = Space(rowpts, colpts)
CreateMesh!(space, rowpts, colpts)
SetDeltas!(space, breadth, length)
fluid = Fluid(rho, mu)


SetTimeStep!(0.8, space)


u_in=1.0 #Lid velocity
v_wall=0.0 #Velocity of fluid at the walls
p_out=1.0
flow=Boundary("D",u_in)
noslip=Boundary("D",v_wall)
zeroflux=Boundary("N",0.0)
pressureatm=Boundary("D",p_out)

setPBoundary!(space, zeroflux, zeroflux, pressureatm, zeroflux)


SolvePressurePoisson!(space, fluid, zeroflux, zeroflux, pressureatm, zeroflux)








setUBoundary!(space,noslip,noslip,flow,noslip)
setVBoundary!(space,noslip,noslip,noslip,noslip)
setPBoundary!(space,noslip,noslip,flow,noslip)


GetStarredVelocities!(space, fluid)


SolvePressurePoisson!(space, fluid, zeroflux, zeroflux, pressureatm, zeroflux)

SolveMomentumEquation!(space,fluid)

























CFL = CFL_number

SetTimeStep!(CFL, space)



setUBoundary!(space, noslip, noslip, flow, noslip)
setVBoundary!(space, noslip, noslip, noslip, noslip)
setPBoundary!(space, zeroflux, zeroflux, pressureatm, zeroflux)
GetStarredVelocities!(space, fluid)

SolvePressurePoisson!(space, water, zeroflux, zeroflux, pressureatm, zeroflux)
SolveMomentumEquation!(space, fluid)
SetCentrePUV!(space)










rows = space.rowpts
cols = space.colpts
u_star = space.u_star
v_star = space.v_star
p = space.p
dx = space.dx
dy = space.dy
dt = space.dt
rho = fluid.rho
factor = 1/(2/dx^2 + 2/dy^2)














# Define the intial error and tolerance for convergence (error > tol at least initially)
error = 1
tol = 1e-3

# Evaluate derivative of starred velocities
u_star1_x = (u_star[2:rows+1, 3:rows+2] - u_star[2:rows+1, 1:cols])/(2*dx)
v_star1_y = (v_star[3:rows+2, 2:rows+1] - v_star[1:rows, 2:cols+1])/(2*dy)

# Continue iterative solution until error becomes smaller than tolerance
i = 0 
while (error>tol)
    i = i + 1
    #print(i)

    # Save current pressure as p_old
    p_old = copy(p)


    # Evaluate second derivative of pressure from p_old
    p2_xy = (p_old[3:rows+2, 2:cols+1] + p_old[1:rows, 2:cols+1])/dy^2 + (p_old[2:rows+1, 3:cols+2] + p_old[2:rows+1, 1:cols])/dx^2

    # Calculate new pressure

    p[2:rows+1, 2:cols+1] = (p2_xy)*factor - (rho*factor/dt)*(u_star1_x+v_star1_y)

    #Find maximum error between old and new pressure matrices
    
    error = maximum(abs.(p-p_old))
    println(error)

    #Apply pressure boundary conditions

    setPBoundary!(space, zeroflux, zeroflux, pressureatm, zeroflux)

    #Escape condition in case solution does not converge after 500 iterations
    if i > 500
        tol = tol*10

    end
end


println("oi")


























































u_star = copy(u)
v_star = copy(v)

    # Calculate the derivatives of u and v using the finite difference scheme

    #u1_y = (u[3:end, 2:end-1] - u[1:end-2, 2:end-1])/(2*dy)
u1_y = (u[3:rows+2, 2:cols+1] - u[1:rows, 2:cols+1])/(2*dy)
u1_x = (u[2:rows+1, 3:cols+2] - u[2:rows+1, 1:cols])/(2*dx)
u2_y = (u[3:rows+2, 2:cols+1] - 2*u[2:rows+1, 2:cols+1] + u[1:rows, 2:cols+1])/(dy^2)
u2_x = (u[2:rows+1, 3:cols+2] - 2*u[2:rows+1, 2:cols+1] + u[2:rows+1, 1:cols])/(dx^2)
v_face = (v[2:rows+1, 2:cols+1] + v[2:rows+1, 1:cols] + v[3:rows+2, 2:cols+1] + v[3:rows+2, 1:cols])/4
u_star[2:rows+1, 2:cols+1] = u[2:rows+1, 2:cols+1] - dt*(u[2:rows+1, 2:cols+1]*u1_x + v_face*u1_y) + (dt*(mu/rho)*(u2_x+u2_y)) .+ dt*S_x

u_star