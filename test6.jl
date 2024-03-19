using Plots
pyplot()


r_ini = 1.0
r_end = 2.0 # range in r-direction
θ_ini = 0.0
θ_end = π/2 # range in θ-direction
n = 10

# Spatial grid for the numerical solution 
r = range(r_ini, stop=r_end, length=n)
θ = range(θ_ini, stop=θ_end, length=n)'

x1 = r .* cos.(θ)
y1 = r .* sin.(θ)
z = ones(n, n).*[1,2,3,4,5,6,7,8,9,10] # Ensure z has the correct dimensions

contourf(x1, y1, z, title="Initial Condition", color=:viridis)