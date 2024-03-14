# Here we list the functions needed to solve the numerical problem:
#-----------------------------------------------------------------------------------------------------------------------------------
# ETBS - Explicit Time Backward Space
#=================#
# Non-vectorized: #
#=================#
# For Dirichlet Boundary Conditions 
# and no comparison with the analytical solution

function scheme(u, un, i, lmbda)
    # ETBS - Explicit Time Backward Space
    u[i+1] = un[i] - lmbda * (un[i] - un[i-1])

    return u
end

function update_soln(nt, nx, lmbda)

    # Numerical solution
    u = copy(uinit)
    u_hist = [copy(u)]

    for tstep in 1:nt
        un = copy(u)
        for i in 2:nx-1
            # We can update the solution in the for loop, or simply call scheme
            #u[i+1] = un[i] - lmbda*(un[i] - un[i-1])
            u = scheme(u, un, i, lmbda)
        end
        u[1] = 0
        u[end] = 0
        push!(u_hist, copy(u))


    end

    return u_hist
end


#=================#
#   Vectorized:   #
#=================#
# For Periodic Boundary Conditons
# and comparison with the analytical solution

function scheme_vectorized(u, x_vec, n, lmbda)
    # Assigning J positions to u we are vectorizing the operations as well
    # that way we do not need a spatial for loop, only temporal

    # ETBS - Explicit Time Backward Space
    u[J] = u[J] - lmbda * (u[J] - u[Jm1])

    #u_exact = f(x_vec .- c*n*ht)
    u_exact = f(mod.(x_vec .- c * (n + 1) * ht, 1.0))


    return u, u_exact
end

function update_soln_vectorized(nt, x_vec, lmbda)

    # Numerical solution
    u = copy(uinit)
    u_hist = [copy(u)]

    # Exact solution
    u_exact = copy(uinit_exact)
    u_hist_exact = [copy(u_exact)]


    for tstep in 1:nt
        u, u_exact = scheme_vectorized(u, x_vec, tstep, lmbda)
        push!(u_hist, copy(u))
        push!(u_hist_exact, copy(u_exact))

    end

    return u_hist, u_hist_exact
end
#-----------------------------------------------------------------------------------------------------------------------------------