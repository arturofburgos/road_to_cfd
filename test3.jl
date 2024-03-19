using Plots, LaTeXStrings
pyplot()

r_ini = 1.0
r_end = 2.0 # range in r-direction
θ_ini = 0.0
θ_end = π/2 # range in θ-direction
n = 100

# Spatial grid for the numerical solution 
r = range(r_ini, stop=r_end, length=n)
θ = range(θ_ini, stop=θ_end, length=n)'

x1 = r * cos.(θ)
y1 = r .* sin.(θ)
z = ones(n, n) # Ensure z has the correct dimensions

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 10


contourf(x1, y1, z, title=L"Plot of $\Gamma_3(x)$",  color=:viridis)