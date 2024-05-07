module FlowJulia

export Boundary, Space, CreateMesh!, SetDeltas!, SetTimeStep!

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
function setUBoundary!(space::Space, left::Boundary, right::Boundary,
    top::Boundary, bottom::Boundary)

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


function setVBoundary!(space::Space, left::Boundary, right::Boundary,
    top::Boundary, bottom::Boundary)

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

function setPBoundary!(space::Space, left::Boundary, right::Boundary,
    top::Boundary, bottom::Boundary)

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


function GetStarredVelocities(space::Space, fluid::Fluid)
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
    u2_y = (u[3:rows+2, 2:cols+1] - 2*u[3:rows+2, 2:cols+1] +
     u[2:rows+1, 1:cols])

    u1_x = (u[2:end-1, 3:end] - u[2:end-1, 1:end-2])/(2*dx)



end