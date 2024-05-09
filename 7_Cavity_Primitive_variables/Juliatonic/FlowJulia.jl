module FlowJulia

using FileIO
export Boundary, Space, CreateMesh!, Fluid, SetCentrePUV!
export SetDeltas!, SetTimeStep!, setPBoundary!, WriteToFile
export GetStarredVelocities!, SolvePressurePoisson!, SolveMomentumEquation!
export setUBoundary!, setVBoundary!, MakeResultDirectory

struct Boundary
    type::String
    value::Float64

    function Boundary(type::String, value::Float64)
        new(type, value)
    end
end

# Create the Space Struct

mutable struct Space
    rowpts::Int
    colpts::Int
    u::Array{Float64,2}
    v::Array{Float64,2}
    u_star::Array{Float64,2}
    v_star::Array{Float64,2}
    u_next::Array{Float64,2}
    v_next::Array{Float64,2}
    u_c::Array{Float64,2}
    v_c::Array{Float64,2}
    p::Array{Float64,2}
    p_c::Array{Float64,2}
    dx::Float64
    dy::Float64
    S_x::Float64
    S_y::Float64
    dt::Float64

    function Space(rowpts, colpts)
        u = zeros(Float64, rowpts + 2, colpts + 2)
        v = zeros(Float64, rowpts + 2, colpts + 2)
        u_star = zeros(Float64, rowpts + 2, colpts + 2)
        v_star = zeros(Float64, rowpts + 2, colpts + 2)
        u_next = zeros(Float64, rowpts + 2, colpts + 2)
        v_next = zeros(Float64, rowpts + 2, colpts + 2)
        u_c = zeros(Float64, rowpts, colpts)
        v_c = zeros(Float64, rowpts, colpts)
        p = zeros(Float64, rowpts + 2, colpts + 2)
        p_c = zeros(Float64, rowpts, colpts)
        dx = 0.0
        dy = 0.0
        S_x = 0.0
        S_y = 0.0
        dt = 0.0

        new(rowpts, colpts, u, v,
            u_star, v_star, u_next, v_next,
            u_c, v_c, p, p_c, dx, dy, S_x, S_y, dt)
    end
end

# Create the inplace modifying functions:

# Create Mesh

function CreateMesh!(space::Space, rowpts::Int, colpts::Int)
    space.rowpts = rowpts
    space.colpts = colpts
    space.u = zeros(Float64, space.rowpts + 2, space.colpts + 2)
    space.v = zeros(Float64, space.rowpts + 2, space.colpts + 2)
    space.u_star = zeros(Float64, space.rowpts + 2, space.colpts + 2)
    space.v_star = zeros(Float64, space.rowpts + 2, space.colpts + 2)
    space.u_next = zeros(Float64, space.rowpts + 2, space.colpts + 2)
    space.v_next = zeros(Float64, space.rowpts + 2, space.colpts + 2)
    space.u_c = zeros(Float64, space.rowpts, space.colpts)
    space.v_c = zeros(Float64, space.rowpts, space.colpts)
    space.p = zeros(Float64, space.rowpts + 2, space.colpts + 2)
    space.p_c = zeros(Float64, space.rowpts, space.colpts)
    SetSourceTerm!(space)
end

# Set dx and dy

function SetDeltas!(space::Space, breadth::Float64, length::Float64)
    space.dx = length / (space.colpts - 1)
    space.dy = breadth / (space.rowpts - 1)
end

# Optional functions for the changing the Space Struct

function SetInitialU!(space::Space, U::Float64)
    space.u .= U .* space.u
end

function SetInitialV!(space::Space, V::Float64)
    space.v .= V .* space.v
end

function SetInitialP!(space::Space, P::Float64)
    space.p .= P .* space.p
end

function SetSourceTerm!(space::Space, S_x::Float64=0.0, S_y::Float64=0.0)
    space.S_x = S_x
    space.S_y = S_y
end



struct Fluid
    rho::Float64
    mu::Float64

    function Fluid(rho::Float64, mu::Float64)
        new(rho, mu)
    end
end


# Set boundary conditions for horizontal velocity
function setUBoundary!(space::Space,
                       left::Boundary,
                       right::Boundary,
                       top::Boundary,
                       bottom::Boundary)

    # Left
    if left.type == "D"
        space.u[:, 1] .= left.value
    elseif left.type == "N"
        space.u[:, 1] .= -left.value * space.dx .+ space.u[:, 2]
    end

    # Right
    if right.type == "D"
        space.u[:, end] .= right.value
    elseif right.type == "N"
        space.u[:, end] .= right.value * space.dx .+ space.u[:, end-1]
    end

    # Top
    if top.type == "D"
        space.u[end, :] .= 2 * top.value .- space.u[end-1, :]
    elseif top.type == "N"
        space.u[end, :] .= -top.value * space.dy .+ space.u[end-1, :]
    end

    # Bottom
    if bottom.type == "D"
        space.u[1, :] .= 2 * bottom.value .- space.u[2, :]
    elseif bottom.type == "N"
        space.u[1, :] .= bottom.value * space.dy .+ space.u[2, :]
    end
end


function setVBoundary!(space::Space,
                       left::Boundary,
                       right::Boundary,
                       top::Boundary,
                       bottom::Boundary)

    # Left
    if left.type == "D"
        space.v[:, 1] .= 2 * left.value .- space.v[:, 2]
    elseif left.type == "N"
        space.v[:, 1] .= -left.value * space.dx .+ space.v[:, 2]
    end

    # Right
    if right.type == "D"
        space.v[:, end] .= 2 * right.value .- space.v[:, end-1]
    elseif right.type == "N"
        space.v[:, end] .= right.value * space.dx .+ space.v[:, end-1]
    end

    # Top
    if top.type == "D"
        space.v[end, :] .= top.value
    elseif top.type == "N"
        space.v[end, :] .= -top.value * space.dy .+ space.v[end-1, :]
    end

    # Bottom
    if bottom.type == "D"
        space.v[1, :] .= bottom.value
    elseif bottom.type == "N"
        space.v[1, :] .= bottom.value * space.dy .+ space.v[2, :]
    end
end

function setPBoundary!(space::Space,
                       left::Boundary,
                       right::Boundary,
                       top::Boundary,
                       bottom::Boundary)

    # Left
    if left.type == "D"
        space.p[:, 1] .= left.value
    elseif left.type == "N"
        space.p[:, 1] .= -left.value * space.dx .+ space.p[:, 2]
    end

    # Right
    if right.type == "D"
        space.p[:, end] .= right.value
    elseif right.type == "N"
        space.p[:, end] .= right.value * space.dx .+ space.p[:, end-1]
    end

    # Top
    if top.type == "D"
        space.p[end, :] .= top.value
    elseif top.type == "N"
        space.p[end, :] .= -top.value * space.dy .+ space.p[end-1, :]
    end

    # Bottom
    if bottom.type == "D"
        space.p[1, :] .= bottom.value
    elseif bottom.type == "N"
        space.p[1, :] .= bottom.value * space.dy .+ space.p[2, :]
    end
end

function SetTimeStep!(CFL::Float64, space::Space)
    dt = CFL / sum([maximum(space.u) / space.dx, maximum(space.v) / space.dy])
    # Escape condition if dt is infinity due to zero velocity initially
    if isinf(dt)
        dt = CFL * (space.dx + space.dy)
    end
    space.dt = dt
end


function GetStarredVelocities!(space::Space, fluid::Fluid)
    # Save struct attributes for improve readability

    rows = space.rowpts
    cols = space.colpts
    u = space.u
    v = space.v
    dx = space.dx
    dy = space.dy
    dt = space.dt
    S_x = space.S_x
    S_y = space.S_y
    rho = fluid.rho
    mu = fluid.mu

    # Copy u and v to the new variables u_star and v_star

    # TODO do I need to perform this copy? Ask Nick

    u_star = copy(u)
    v_star = copy(v)

    # Calculate the derivatives of u and v using the finite difference scheme

    #u1_y = (u[3:end, 2:end-1] - u[1:end-2, 2:end-1])/(2*dy)
    u1_y = (u[3:rows+2, 2:cols+1] - u[1:rows, 2:cols+1])/(2*dy)
    u1_x = (u[2:rows+1, 3:cols+2] - u[2:rows+1, 1:cols])/(2*dx)
    u2_y = (u[3:rows+2, 2:cols+1] - 2*u[2:rows+1, 2:cols+1] +
     u[1:rows, 2:cols+1])
    u2_x = (u[3:rows+2, 2:cols+1] - 2*u[2:rows+1, 2:cols+1] +
     u[1:rows, 2:cols+1])
    v_face = (v[2:rows+1, 2:cols+1] + v[2:rows+1, 1:cols] +
     v[3:rows+2, 2:cols+1] + v[3:rows+2, 1:cols])/4
    u_star[2:rows+1, 2:cols+1] = u[2:rows+1, 2:cols+1] 
     - dt*(u[2:rows+1, 2:cols+1]*u1_x + v_face*u1_y) 
     + dt*(mu/rho)*(u2_x+u2_y) .+ dt*S_x


    v1_y = (v[3:rows+2, 2:cols+1] - v[1:rows, 2:cols+1])/(2*dy)
    v1_x = (v[2:rows+1, 3:cols+2] - v[2:rows+1, 1:cols])/(2*dx)
    v2_y = (v[3:rows+2, 2:cols+1] - 2*v[2:rows+1, 2:cols+1] +
     v[1:rows, 2:cols+1])
    v2_x = (v[3:rows+2, 2:cols+1] - 2*v[2:rows+1, 2:cols+1] +
     v[1:rows, 2:cols+1])
    u_face = (u[2:rows+1, 2:cols+1] + u[2:rows+1, 3:cols+2] +
     u[1:rows, 2:cols+1] + u[1:rows, 3:cols+2])/4
    v_star[2:rows+1, 2:cols+1] = v[2:rows+1, 2:cols+1] 
     - dt*(v[2:rows+1, 2:cols+1]*v1_y + u_face*v1_x) 
     + dt*(mu/rho)*(v2_x+v2_y) .+ dt*S_y

    # Save the calculated starred velocities to the space struct

    space.u_star = copy(u_star)
    space.v_star = copy(v_star)

end

# The second function is used to iteratively solve the pressure Possion equation from the starred velocities 
# to calculate pressure at t+delta_t


function SolvePressurePoisson!(space::Space,
                              fluid::Fluid,
                              left::Boundary,
                              right::Boundary,
                              top::Boundary,
                              bottom::Boundary)

    # Save struct attributes for improve readability

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
    v_star1_y = (v_star[3:rows+2, 2:rows+1] - u_star[1:rows, 2:cols+1])/(2*dy)
 
    # Continue iterative solution until error becomes smaller than tolerance
    i = 0 
    while (error>tol)
        i = i + 1
        #print(i)

        # Save current pressure as p_old
        p_old = copy(p)


        # Evaluate second derivative of pressure from p_old
        p2_xy = (p_old[3:rows+2, 2:cols+1] + p_old[1:rows, 2:cols+1])/dy^2 +
         (p_old[2:rows+1, 3:cols+2] + p_old[2:rows+1, 1:cols])/dx^2

        # Calculate new pressure

        p[2:rows+1, 2:cols+1] = (p2_xy)*factor - (rho*factor/dt)*(u_star1_x+v_star1_y)

        #Find maximum error between old and new pressure matrices
        
        error = maximum(abs.(p-p_old))

        #Apply pressure boundary conditions

        setPBoundary!(space, left, right, top, bottom)

        #Escape condition in case solution does not converge after 500 iterations
        if i > 500
            tol = tol*10

        end
    end
end


function SolveMomentumEquation!(space::Space, fluid::Fluid)
    # Save struct attributes for improve readability

    rows = space.rowpts
    cols = space.colpts
    u_star = space.u_star
    v_star = space.v_star
    p = space.p
    dx = space.dx
    dy = space.dy
    dt = space.dt
    rho = fluid.rho
    u = space.u
    v = space.v

    #Evaluate first derivative of pressure in x direction
    p1_x = (p[2:rows+1, 3:rows+2] - p[2:rows+1, 1:cols])/(2*dx)

    # Calculate u at the next SetTimeStep
    u[2:rows+1, 2:cols+1]= u_star[2:rows+1, 2:cols+1]-(dt/rho)*p1_x

    #Evaluate first derivative of pressure in y direction
    p1_y = (p[3:rows+2, 2:rows+1] - p[1:rows, 2:cols+1])/(2*dy)

    #Calculate v at next timestep
    v[2:rows+1, 2:cols+1]= v_star[2:rows+1,2:cols+1]-(dt/rho)*p1_y

end


function SetCentrePUV!(space::Space)
    space.p_c = space.p[2:end-1, 2:end]
    space.u_c = space.u[2:end-1, 2:end]
    space.v_c = space.v[2:end-1, 2:end]
end

function MakeResultDirectory(;wipe::Bool=true)

    # Get path to the Result directory
    cwdir = pwd()
    inter_path = joinpath(cwdir, "7_Cavity_Primitive_variables", "Juliatonic")
    dir_path = joinpath(inter_path, "Result")
    

    # If directory does not exists, make iterations
    
    if !isdir(dir_path)
        mkdir(dir_path)
    else
        # If wipe is true, remove files present  in the directory
        if wipe
            cd(dir_path)
            filelist = readdir()
            for file in filelist
                rm(file)
            end
        end
    end
    cd(cwdir)

end

function WriteToFile(space::Space, iteration, interval)
    if (1/2 * iteration % interval == 0)
        cwdir = pwd()
        inter_path = joinpath(cwdir, "7_Cavity_Primitive_variables", "Juliatonic")
        dir_path = joinpath(inter_path, "Result") 
        filename = "PUV$(iteration).txt"
        path = joinpath(dir_path, filename)
        open(path, "w") do f
            for i in 1:space.rowpts
                for j in 1:space.colpts
                    println(f, "$(space.p_c[i, j])\t$(space.u_c[i, j])\t$(space.v_c[i, j])")
                end
            end
        end
    end
end


end