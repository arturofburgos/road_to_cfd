using Plots

pythonplot()
s = range(1, 2, length=100)
phi = range(0, 2π, length=100)

x = [s * cos(p) for s in s, p in phi]
y = [s * sin(p) for s in s, p in phi]
z = ones(size(x)) 

contourf(x, y, z)
display(gcf())