using LaTeXStrings, Plots
#pyplot()

f(x, y) = (3x + y^2) * abs(sin(x) + cos(y))

x = range(0, 5, length=100)
y = range(0, 3, length=50)
z = @. f(x', y)

contour(x, y, z, levels=10, color=:viridis, clabels=true, cbar=false, lw=1)
title!(L"Plot of $(3x + y^2)|\sin(x) + \cos(y)|$")
xlabel!(L"x")
ylabel!(L"y")

contour(x, y, z, color=[:black])

contour(x, y, z, fill=true)


contourf(x, y, z, levels=20, color=:seaborn_rocket_gradient)
title!(L"(3x + y^2)|\sin(x) + \cos(y)|")
xlabel!(L"x")
ylabel!(L"y")

heatmap(x, y, z, levels=20, color=:seaborn_rocket_gradient)
heatmap(x, y, z, color=:seaborn_rocket_gradient)