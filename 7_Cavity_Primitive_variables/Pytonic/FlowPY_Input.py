import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from FlowPY import *

class Space:
    def __init__(self):
        pass
    
    def CreateMesh(self,rowpts,colpts):
        #Domain gridpoints
        self.rowpts=rowpts
        self.colpts=colpts        #Velocity matrices
        self.u=np.zeros((self.rowpts+2,self.colpts+2))
        self.v=np.zeros((self.rowpts+2,self.colpts+2))
        self.u_star=np.zeros((self.rowpts+2,self.colpts+2))
        self.v_star=np.zeros((self.rowpts+2,self.colpts+2))
        self.u_next=np.zeros((self.rowpts+2,self.colpts+2))
        self.v_next=np.zeros((self.rowpts+2,self.colpts+2))
        self.u_c=np.zeros((self.rowpts,self.colpts))
        self.v_c=np.zeros((self.rowpts,self.colpts))
        #Pressure matrices
        self.p=np.zeros((self.rowpts+2,self.colpts+2))
        self.p_c=np.zeros((self.rowpts,self.colpts))
        #Set default source term
        self.SetSourceTerm()        
        
    def SetDeltas(self,breadth,length):
        self.dx=length/(self.colpts-1)
        self.dy=breadth/(self.rowpts-1)    
        
    def SetInitialU(self,U):
        self.u=U*self.u
        
    def SetInitialV(self,V):
        self.v=V*self.v
        
    def SetInitialP(self,P):
        self.p=P*self.p    
        
    def SetSourceTerm(self,S_x=0,S_y=0):
        self.S_x=S_x
        self.S_y=S_y
#### SPATIAL AND TEMPORAL INPUTS
length=4 #Length of computational domain in the x-direction
breadth=4 #Breadth of computational domain in the y-direction
colpts=10 #Number of grid points in the x-direction #KEEP ODD
rowpts = 10 #Number of grid points in the y-direction #KEEP ODD#Create an object of the class Space called cavity
cavity=Space()
cavity.CreateMesh(rowpts,colpts)
cavity.SetDeltas(breadth,length)

class Fluid:
    def __init__(self,rho,mu):
        self.SetFluidProperties(rho,mu)
    
    def SetFluidProperties(self,rho,mu):
        self.rho=rho
        self.mu=mu
#### FLUID PROPERTIES
rho=1 #Density of fluid
mu=0.01 #Viscosity of fluid#Create an object of the class Fluid called water
water=Fluid(rho,mu)

def SetTimeStep(CFL,space,fluid):
    with np.errstate(divide='ignore'):
        dt=CFL/np.sum([np.amax(space.u)/space.dx,\
                           np.amax(space.v)/space.dy])
    #Escape condition if dt is infinity due to zero velocity initially
    if np.isinf(dt):
        dt=CFL*(space.dx+space.dy)
    space.dt=dt



SetTimeStep(0.8, cavity, water)


def GetStarredVelocities(space,fluid):
    #Save object attributes as local variable with explicity typing for improved readability
    rows=int(space.rowpts)
    cols=int(space.colpts)
    u=space.u.astype(float,copy=False)
    v=space.v.astype(float,copy=False)
    dx=float(space.dx)
    dy=float(space.dy)
    dt=float(space.dt)
    S_x=float(space.S_x)
    S_y=float(space.S_y)
    rho=float(fluid.rho)
    mu=float(fluid.mu)

    
    
    #Copy u and v to new variables u_star and v_star
    u_star=u.copy()
    v_star=v.copy()
    
    #Calculate derivatives of u and v using the finite difference scheme. 
    #Numpy vectorization saves us from using slower for loops to go over each element in the u and v matrices
    u1_y=(u[2:,1:cols+1]-u[0:rows,1:cols+1])/(2*dy)
    u1_x=(u[1:rows+1,2:]-u[1:rows+1,0:cols])/(2*dx)
    u2_y=(u[2:,1:cols+1]-2*u[1:rows+1,1:cols+1]+u[0:rows,1:cols+1])/(dy**2)
    u2_x=(u[1:rows+1,2:]-2*u[1:rows+1,1:cols+1]+u[1:rows+1,0:cols])/(dx**2)
    v_face=(v[1:rows+1,1:cols+1]+v[1:rows+1,0:cols]+v[2:,1:cols+1]+v[2:,0:cols])/4
    u_star[1:rows+1,1:cols+1]=u[1:rows+1,1:cols+1]-dt*(u[1:rows+1,1:cols+1]*u1_x+v_face*u1_y)+(dt*(mu/rho)*(u2_x+u2_y))+(dt*S_x)   

    v1_y=(v[2:,1:cols+1]-v[0:rows,1:cols+1])/(2*dy)
    v1_x=(v[1:rows+1,2:]-v[1:rows+1,0:cols])/(2*dx)
    v2_y=(v[2:,1:cols+1]-2*v[1:rows+1,1:cols+1]+v[0:rows,1:cols+1])/(dy**2)
    v2_x=(v[1:rows+1,2:]-2*v[1:rows+1,1:cols+1]+v[1:rows+1,0:cols])/(dx**2)
    u_face=(u[1:rows+1,1:cols+1]+u[1:rows+1,2:]+u[0:rows,1:cols+1]+u[0:rows,2:])/4
    v_star[1:rows+1,1:cols+1]=v[1:rows+1,1:cols+1]-dt*(u_face*v1_x+v[1:rows+1,1:cols+1]*v1_y)+(dt*(mu/rho)*(v2_x+v2_y))+(dt*S_y)
    
    #Save the calculated starred velocities to the space object 
    space.u_star=u_star.copy()
    space.v_star=v_star.copy()



def SolvePressurePoisson(space,fluid,left,right,top,bottom):
    #Save object attributes as local variable with explicity typing for improved readability
    rows=int(space.rowpts)
    cols=int(space.colpts)
    u_star=space.u_star.astype(float,copy=False)
    v_star=space.v_star.astype(float,copy=False)
    p=space.p.astype(float,copy=False)
    dx=float(space.dx)
    dy=float(space.dy)
    dt=float(space.dt)
    rho=float(fluid.rho)
    factor=1/(2/dx**2+2/dy**2)
    
    #Define initial error and tolerance for convergence (error > tol necessary initially)
    error=1
    tol=1e-3

    #Evaluate derivative of starred velocities
    ustar1_x=(u_star[1:rows+1,2:]-u_star[1:rows+1,0:cols])/(2*dx)
    vstar1_y=(v_star[2:,1:cols+1]-v_star[0:rows,1:cols+1])/(2*dy)

    #Continue iterative solution until error becomes smaller than tolerance
    i=0
    while(error>tol):
        i+=1
        
        #Save current pressure as p_old
        p_old=p.astype(float,copy=True)
        
        #Evaluate second derivative of pressure from p_old
        p2_xy=(p_old[2:,1:cols+1]+p_old[0:rows,1:cols+1])/dy**2+(p_old[1:rows+1,2:]+p_old[1:rows+1,0:cols])/dx**2
        
        #Calculate new pressure 
        p[1:rows+1,1:cols+1]=(p2_xy)*factor-(rho*factor/dt)*(ustar1_x+vstar1_y)
        
        #Find maximum error between old and new pressure matrices
        error=np.amax(abs(p-p_old))
        
        #Apply pressure boundary conditions
        SetPBoundary(space,left,right,top,bottom)
        
        #Escape condition in case solution does not converge after 500 iterations
        if(i>500):
            tol*=10

GetStarredVelocities(cavity, water)

class Boundary:
    def __init__(self,boundary_type,boundary_value):
        self.DefineBoundary(boundary_type,boundary_value)
        
    def DefineBoundary(self,boundary_type,boundary_value):
        self.type=boundary_type
        self.value=boundary_value

#### BOUNDARY SPECIFICATIONS
u_in=1 #Lid velocity
v_wall=0 #Velocity of fluid at the walls
p_out=1 #Gauge pressure at the boundaries#Create objects of the class Boundary having either Dirichlet ("D") or Neumann ("N") type boundaries
flow=Boundary("D",u_in)
noslip=Boundary("D",v_wall)
zeroflux=Boundary("N",0)
pressureatm=Boundary("D",p_out)

def SetPBoundary(space,left,right,top,bottom):
    if(left.type=="D"):
        space.p[:,0]=left.value
    elif(left.type=="N"):
        space.p[:,0]=-left.value*space.dx+space.p[:,1]
    
    if(right.type=="D"):
        space.p[:,-1]=right.value # Error maybe ? it was p[1,-1]
    elif(right.type=="N"):
        space.p[:,-1]=right.value*space.dx+space.p[:,-2]
        
    if(top.type=="D"):
        space.p[-1,:]=top.value
    elif(top.type=="N"):
        space.p[-1,:]=-top.value*space.dy+space.p[-2,:]
     
    if(bottom.type=="D"):
        space.p[0,:]=bottom.value
    elif(bottom.type=="N"):
        space.p[0,:]=bottom.value*space.dy+space.p[1,:]

SetPBoundary(cavity, zeroflux, zeroflux, pressureatm, zeroflux)

SolvePressurePoisson(cavity,water,zeroflux,zeroflux, pressureatm, zeroflux)    
    


def SetUBoundary(space,left,right,top,bottom):
    if(left.type=="D"):
        space.u[:,0]=left.value
    elif(left.type=="N"):
        space.u[:,0]=-left.value*space.dx+space.u[:,1]
    
    if(right.type=="D"):
        space.u[:,-1]=right.value
    elif(right.type=="N"):
        space.u[:,-1]=right.value*space.dx+space.u[:,-2]
        
    if(top.type=="D"):
        space.u[-1,:]=2*top.value-space.u[-2,:]
    elif(top.type=="N"):
        space.u[-1,:]=-top.value*space.dy+space.u[-2,:]
     
    if(bottom.type=="D"):
        space.u[0,:]=2*bottom.value-space.u[1,:]
    elif(bottom.type=="N"):
        space.u[0,:]=bottom.value*space.dy+space.u[1,:]

#Set boundary conditions for vertical velocity
def SetVBoundary(space,left,right,top,bottom):
    if(left.type=="D"):
        space.v[:,0]=2*left.value-space.v[:,1]
    elif(left.type=="N"):
        space.v[:,0]=-left.value*space.dx+space.v[:,1]
    
    if(right.type=="D"):
        space.v[:,-1]=2*right.value-space.v[:,-2]
    elif(right.type=="N"):
        space.v[:,-1]=right.value*space.dx+space.v[:,-2]
        
    if(top.type=="D"):
        space.v[-1,:]=top.value
    elif(top.type=="N"):
        space.v[-1,:]=-top.value*space.dy+space.v[-2,:]
     
    if(bottom.type=="D"):
        space.v[0,:]=bottom.value
    elif(bottom.type=="N"):
        space.v[0,:]=bottom.value*space.dy+space.v[1,:]


def SolveMomentumEquation(space,fluid):
    #Save object attributes as local variable with explicity typing for improved readability
    rows=int(space.rowpts)
    cols=int(space.colpts)
    u_star=space.u_star.astype(float,copy=False)
    v_star=space.v_star.astype(float,copy=False)
    p=space.p.astype(float,copy=False)
    dx=float(space.dx)
    dy=float(space.dy)
    dt=float(space.dt)
    rho=float(fluid.rho)
    u=space.u.astype(float,copy=False)
    v=space.v.astype(float,copy=False)

    #Evaluate first derivative of pressure in x direction
    p1_x=(p[1:rows+1,2:]-p[1:rows+1,0:cols])/(2*dx)
    #Calculate u at next timestep
    u[1:rows+1,1:cols+1]=u_star[1:rows+1,1:cols+1]-(dt/rho)*p1_x

    #Evaluate first derivative of pressure in y direction
    p1_y=(p[2:,1:cols+1]-p[0:rows,1:cols+1])/(2*dy)
    #Calculate v at next timestep
    v[1:rows+1,1:cols+1]=v_star[1:rows+1,1:cols+1]-(dt/rho)*p1_y


def SetCentrePUV(space):
    space.p_c=space.p[1:-1,1:-1]
    space.u_c=space.u[1:-1,1:-1]
    space.v_c=space.v[1:-1,1:-1]


SetTimeStep(0.8, cavity, water)


SetUBoundary(cavity,noslip,noslip,flow,noslip)

SetVBoundary(cavity,noslip,noslip,noslip,noslip)
SetPBoundary(cavity,zeroflux,zeroflux,pressureatm,zeroflux)


GetStarredVelocities(cavity,water)

SolvePressurePoisson(cavity,water,zeroflux,zeroflux, pressureatm, zeroflux)    


SolveMomentumEquation(cavity,water)
SetCentrePUV(cavity)





dt = cavity.dt
dx = cavity.dx 

dy = cavity.dy 

rows = cavity.rowpts
cols = cavity.colpts
u = cavity.u

u_star = u.copy()

u1_y=(u[2:,1:cols+1]-u[0:rows,1:cols+1])/(2*dy)
u1_x=(u[1:rows+1,2:]-u[1:rows+1,0:cols])/(2*dx)
u2_y=(u[2:,1:cols+1]-2*u[1:rows+1,1:cols+1]+u[0:rows,1:cols+1])/(dy**2)
u2_x=(u[1:rows+1,2:]-2*u[1:rows+1,1:cols+1]+u[1:rows+1,0:cols])/(dx**2)
#v_face=(v[1:rows+1,1:cols+1]+v[1:rows+1,0:cols]+v[2:,1:cols+1]+v[2:,0:cols])/4
u_star[1:rows+1,1:cols+1]=u[1:rows+1,1:cols+1]-dt*(u[1:rows+1,1:cols+1]*u1_x+v_face*u1_y)+(dt*(mu/rho)*(u2_x+u2_y))+(dt*S_x)   






rows = cavity.rowpts
cols = cavity.colpts
u_star = cavity.u_star
v_star = cavity.v_star
p = cavity.p
dx = cavity.dx
dy = cavity.dy
dt = cavity.dt
rho = water.rho
factor = 1/(2/dx**2 + 2/dy**2)

#Define initial error and tolerance for convergence (error > tol necessary initially)
error=1
tol=1e-3

#Evaluate derivative of starred velocities
ustar1_x=(u_star[1:rows+1,2:]-u_star[1:rows+1,0:cols])/(2*dx)
vstar1_y=(v_star[2:,1:cols+1]-v_star[0:rows,1:cols+1])/(2*dy)
#Continue iterative solution until error becomes smaller than tolerance
i=0
while(error>tol):
    i+=1
    
    #Save current pressure as p_old
    p_old=p.astype(float,copy=True)
    
    #Evaluate second derivative of pressure from p_old
    p2_xy=(p_old[2:,1:cols+1]+p_old[0:rows,1:cols+1])/dy**2+(p_old[1:rows+1,2:]+p_old[1:rows+1,0:cols])/dx**2
    
    #Calculate new pressure 
    p[1:rows+1,1:cols+1]=(p2_xy)*factor-(rho*factor/dt)*(ustar1_x+vstar1_y)
    
    #Find maximum error between old and new pressure matrices
    error=np.amax(abs(p-p_old))
    
    print(error)
    print("\n")
    #Apply pressure boundary conditions
    SetPBoundary(cavity,zeroflux, zeroflux, pressureatm, zeroflux)
    
    #Escape condition in case solution does not converge after 500 iterations
    if(i>500):
        tol*=10


print("oi")




#### SIMULATION PARAMETERS
time=5 #Simulation time
CFL_number=0.5 #Reduce this if solution diverges
file_flag=1 #Keep 1 to print results to file
interval=100 #Record values in file per interval number of iterations
plot_flag=1 #Keep 1 to plot results at the end

#### RUN SIMULATION# Print general simulation information
print("######## Beginning FlowPy Simulation ########")
print("#############################################")
print("# Simulation time: {0:.2f}".format(time))
print("# Mesh: {0} x {1}".format(colpts,rowpts))
print("# Re/u: {0:.2f}\tRe/v:{1:.2f}".format(rho*length/mu,rho*breadth/mu))
print("# Save outputs to text file: {0}".format(bool(file_flag)))## Initialization
# Make directory to store results
MakeResultDirectory(wipe=True)# Initialize counters
t=0
i=0## Run
while(t<time):
    #Print time left
    sys.stdout.write("\rSimulation time left: {0:.2f}".format(time-t))
    sys.stdout.flush()    #Set the time-step
    SetTimeStep(CFL_number,cavity,water)
    timestep=cavity.dt
    
    #Set boundary conditions
    SetUBoundary(cavity,noslip,noslip,flow,noslip)
    SetVBoundary(cavity,noslip,noslip,noslip,noslip)
    SetPBoundary(cavity,zeroflux,zeroflux,pressureatm,zeroflux)
    #Calculate starred velocities
    GetStarredVelocities(cavity,water)
    
    #Solve the pressure Poisson equation
    SolvePressurePoisson(cavity,water,zeroflux,zeroflux, pressureatm, zeroflux)    
    #SSolve the momentum equation
    SolveMomentumEquation(cavity,water)    
    #Save variables and write to file
    SetCentrePUV(cavity)
    if(file_flag==1):
        WriteToFile(cavity,i,interval)    #Advance time-step and counter
    t+=timestep
    i+=1


#######################SET ARRAYS FOR PLOTTING#################################
x=np.linspace(0,length,colpts)
y=np.linspace(0,breadth,rowpts)
[X,Y]=np.meshgrid(x,y)

u=cavity.u
v=cavity.v
p=cavity.p
u_c=cavity.u_c
v_c=cavity.v_c
p_c=cavity.p_c

#Ghia et al. Cavity test benchmark
y_g=[0,0.0547,0.0625,0.0703,0.1016,0.1719,0.2813,0.4531,0.5,0.6172,0.7344,0.8516,0.9531,0.9609,0.9688,0.9766]
u_g=[0,-0.08186,-0.09266,-0.10338,-0.14612,-0.24299,-0.32726,-0.17119,-0.11477,0.02135,0.16256,0.29093,0.55892,0.61756,0.68439,0.75837]

x_g=[0,0.0625,0.0703,0.0781,0.0983,0.1563,0.2266,0.2344,0.5,0.8047,0.8594,0.9063,0.9453,0.9531,0.9609,0.9688]
v_g=[0,0.1836,0.19713,0.20920,0.22965,0.28124,0.30203,0.30174,0.05186,-0.38598,-0.44993,-0.23827,-0.22847,-0.19254,-0.15663,-0.12146]

y_g=[breadth*y_g[i] for i in range(len(y_g))]
x_g=[length*x_g[i] for i in range(len(x_g))]


######################EXTRA PLOTTING CODE BELOW################################
if(plot_flag==1):
    plt.figure(figsize=(20,20))
    plt.contourf(X,Y,p_c,cmap=cm.viridis)
    plt.colorbar()
    plt.quiver(X,Y,u_c,v_c)
    plt.title("Velocity and Pressure Plot")
    
    plt.figure(figsize=(20,20))
    plt.plot(y,u_c[:,int(np.ceil(colpts/2))],"darkblue")
    plt.plot(y_g,u_g,"rx")
    plt.xlabel("Vertical distance along center")
    plt.ylabel("Horizontal velocity")
    plt.title("Benchmark plot 1")
    
    plt.figure(figsize=(20,20))
    plt.plot(x,v_c[int(np.ceil(rowpts/2)),:],"darkblue")
    plt.plot(x_g,v_g,"rx")
    plt.xlabel("Horizontal distance along center")
    plt.ylabel("Vertical velocity")
    plt.title("Benchmark plot 2")
    plt.show()
    
#########################END OF FILE###########################################