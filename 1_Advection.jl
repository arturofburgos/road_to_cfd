using Plots, LaTeXStrings, Measures # Possible to use PyPlot too.
gr()

# Plot defaults
plot_font = "Computer Modern"
default(fontfamily = plot_font, linewidth = 2,
framestyle = :box, grid = false, gridlinewidth = 1,
gridalpha = 0.1, minorgrid = false, label = false, 
ylim = (-0.1,2.1))


c = 1.0
T = 1.0/c

nx = 65 # n + 1 = 65 -> n = 64

x = range(0,1,nx)
hx = x[2] - x[1]
xx = range(0,1,101)
hxx = xx[2] - xx[1]

function square(x)
    u = zeros(size(x))
    u[findall(x -> x >0.4 && x<0.6, x)] .= 2
    return u
end

f = square

uinit = f(x)
uinit_exact = f(xx)

plot(x, uinit, lw = 3, title = "Initial Condition")

lmbda = 0.95
ht = hx*lmbda/c
nt = round(T/ht)


println("     T = $T")
println("tsteps = $nt")
println("    hx = $hx")
println("    ht = $ht")
println("lambda = $lmbda")


J = collect(1:nx)
Jm1 = circshift(J, 1)
Jp1 = circshift(J, -1)

println("J: ", J)
println("Jm1: ", Jm1)
println("Jp1: ", Jp1)

function scheme(u, x_vec, n)
    # Assigning J positions to u we are vectorizing the operations as well
    # that way we do not need a spatial for loop, only temporal

    # ETBS - Explicit Time Backward Space
    u[J] = u[J] - lmbda*(u[J] - u[Jm1])

    #u_exact = f(x_vec .- c*n*ht)
    u_exact = f(mod.(x_vec .- c* (n + 1)* ht, 1.0))
    

    return u, u_exact
end

function update_soln(nt, x_vec)
    
    # Numerical solution
    u = copy(uinit)
    u_hist = [copy(u)]

    # Exact solution
    u_exact = copy(uinit_exact)
    u_hist_exact = [copy(u_exact)]


    for tstep in 1:nt
        u, u_exact = scheme(u, x_vec, tstep)
        push!(u_hist,copy(u))
        push!(u_hist_exact,copy(u_exact))
        
    end

    return u_hist, u_hist_exact
end


function update_plot(n, u_hist, x)
    plot(x, u_hist[n])
    plot!(xx, u_hist_exact[n])
    
end


u_hist, u_hist_exact = update_soln(nt,xx)

anim = @animate for n in 1:length(u_hist)
    update_plot(n, u_hist, x)
end

gif(anim, "ha.gif", fps =15)

