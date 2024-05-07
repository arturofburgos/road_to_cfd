include("FlowJulia.jl")

using .FlowJulia


length = 4.0  # Length of computational domain in the x-direction
breadth = 4.0  # Breadth of computational domain in the y-direction
colpts = 10  # Number of grid points in the x-direction #KEEP ODD
rowpts = 10  # Number of grid points in the y-direction #KEEP ODD

space = Space(rowpts, colpts)
CreateMesh!(space, rowpts, colpts)
SetDeltas!(space, breadth, length)

SetTimeStep!(0.8, space)

space.u