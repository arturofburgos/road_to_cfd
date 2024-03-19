# Here we list the functions needed to solve the numerical problem:
#-----------------------------------------------------------------------------------------------------------------------------------
# Heat Equation 2D:
#=================#
#   Vectorized:   #
#=================#
# For Dirichlet Boundary Conditions 
# and no comparison with the analytical solution

function apply_BC(u)

    u[:, end] .= u_top
    u[:, 1] .= u_bottom
    u[1, :] .= u_left  
    u[end, :] .= u_right

    return u
end

function laplace(u, α, hx, hy)
    rhs = zeros(size(u))

    rhs = α * (((u[Im1, :] - 2 * u[I, :] + u[Ip1, :]) / (hx^2)) +
               ((u[:, Jm1] - 2 * u[:, J] + u[:, Jp1]) / (hy^2)))

    return rhs
end

function update_soln(nt)

    # Numerical solution
    u = copy(u0)
    u_hist = [copy(u)]

    # # Exact solution
    # u_exact = copy(u0)
    # u_hist_exact = [copy(u_exact)]

    for tstep in 1:nt
        #time = tstep*nt # implement this is file hdf5
        # Compute the RHS (2D Laplace)
        rhs = laplace(u, a, hx, hy)

        # Update solution
        u .= u + rhs * ht
        u = apply_BC(u)


        push!(u_hist, copy(u))


    end

    return u_hist

end
#-----------------------------------------------------------------------------------------------------------------------------------