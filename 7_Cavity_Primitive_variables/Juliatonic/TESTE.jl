# Adaptation of the Python code:

mutable struct Boundary
    boundary_type::String
    boundary_value::Float64

    function Boundary(boundary_type::String, boundary_value::Float64)
        new(boundary_type, boundary_value)
    end
end

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

    function Space(rowpts::Int, colpts::Int)
        new(rowpts, colpts,
            zeros(Float64, rowpts+2, colpts+2),
            zeros(Float64, rowpts+2, colpts+2),
            zeros(Float64, rowpts+2, colpts+2),
            zeros(Float64, rowpts+2, colpts+2),
            zeros(Float64, rowpts+2, colpts+2),
            zeros(Float64, rowpts+2, colpts+2),
            zeros(Float64, rowpts, colpts),
            zeros(Float64, rowpts, colpts),
            zeros(Float64, rowpts+2, colpts+2),
            zeros(Float64, rowpts, colpts),
            0.0, 0.0, 0.0, 0.0)
    end

end

function CreateMesh(space::Space, rowpts::Int, colpts::Int)
    space.rowpts = rowpts
    space.colpts = colpts
    space.u = zeros(Float64, rowpts+2, colpts+2)
    space.v = zeros(Float64, rowpts+2, colpts+2)
    space.u_star = zeros(Float64, rowpts+2, colpts+2)
    space.v_star = zeros(Float64, rowpts+2, colpts+2)
    space.u_next = zeros(Float64, rowpts+2, colpts+2)
    space.v_next = zeros(Float64, rowpts+2, colpts+2)
    space.u_c = zeros(Float64, rowpts, colpts)
    space.v_c = zeros(Float64, rowpts, colpts)
    space.p = zeros(Float64, rowpts+2, colpts+2)
    space.p_c = zeros(Float64, rowpts, colpts)
    SetSourceTerm!(space)
    return space
end

function SetDeltas!(space::Space, breadth::Float64, length::Float64)
    space.dx = length / (space.colpts - 1)
    space.dy = breadth / (space.rowpts - 1)
end

function SetInitialU(space::Space, U::Float64)
    space.u .= U .* space.u
end

function SetInitialV(space::Space, V::Float64)
    space.v .= V .* space.v
end

function SetInitialP(space::Space, P::Float64)
    space.p .= P .* space.p
end

function SetSourceTerm!(space::Space; S_x::Float64=0.0, S_y::Float64=0.0)
    space.S_x = S_x
    space.S_y = S_y
end




