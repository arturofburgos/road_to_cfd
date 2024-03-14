# Here we list the functions needed to solve the numerical problem:
#-----------------------------------------------------------------------------------------------------------------------------------
# ETCS - Explicit Time Central Space
#=================#
# Non-vectorized: #
#=================#
# For Dirichlet Boundary Conditions 
# and no comparison with the analytical solution

function scheme(u, un, i, sigma)
    # ETCS - Explicit Time Central Space
    u[i] = un[i] + sigma * (un[i+1] - 2 * un[i] + un[i-1])

    return u
end

function update_soln(nt, nx, sigma)

    # Numerical solution
    u = copy(uinit)
    u_hist = [copy(u)]

    for tstep in 1:nt
        un = copy(u)

        for i in 2:nx-2
            # We can update the solution in the for loop, or simply call scheme
            #u[i] = un[i] + sigma * (un[i+1] - 2*un[i] + un[i-1])
            u = scheme(u, un, i, sigma)
        end

        # Enforcing BCs
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

function scheme_vectorized1(u, sigma)
    # Assigning J positions to u we are vectorizing the operations as well
    # that way we do not need a spatial for loop, only temporal

    # ETCS - Explicit Time Central Space
    u[J] = u[J] + sigma * (u[Jp1] - 2 * u[J] + u[Jm1])

    return u
end

function update_soln_vectorized1(nt, sigma)

    # Numerical solution
    u = copy(uinit)
    u_hist = [copy(u)]

    for tstep in 1:nt
        u = scheme_vectorized1(u, sigma)

        # Enforce dirichlet boundary conditions at the free ends
        # If you comment that next two lines, you have "periodic BCs"
        u[1] = 0
        u[end] = 0
        push!(u_hist, copy(u))

    end

    return u_hist
end


function scheme_vectorized2(u, un, sigma)
    # ETCS - Explicit Time Central Space
    u[2:end-1] = un[2:end-1] + sigma * (un[3:end] - 2 * un[2:end-1] + un[1:end-2])

    return u
end

function update_soln_vectorized2(nt, nx, sigma)

    # Numerical solution
    u = copy(uinit)
    u_hist = [copy(u)]

    for tstep in 1:nt
        un = copy(u)

        u = scheme_vectorized2(u, un, sigma)

        u[1] = 0
        u[end] = 0
        push!(u_hist, copy(u))

    end

    return u_hist
end
#-----------------------------------------------------------------------------------------------------------------------------------