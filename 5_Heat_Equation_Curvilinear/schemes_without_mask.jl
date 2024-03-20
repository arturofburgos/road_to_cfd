# Here we list the functions needed to solve the numerical problem:
#-----------------------------------------------------------------------------------------------------------------------------------
# Heat Equation 2D:
#=================#
#   Vectorized:   #
#=================#
# For Dirichlet Boundary Conditions 
# and no comparison with the analytical solution

function apply_BC(u)

    u[1, :] .= u_left
    u[end, :] .= u_right
    u[:, end] .= u_top
    u[:, 1] .= u_bottom

    return u
end

function laplace(u, α, r, hr, hθ)
    rhs = zeros(size(u))

    rhs = α * (((u[Im1, :] - 2 * u[I, :] + u[Ip1, :]) / (hr^2)) +
               (1 ./ r) .* ((u[Ip1, :] - u[Im1, :]) / (2 * hr)) +
               (1 ./ r .^ 2) .* ((u[:, Jm1] - 2 * u[:, J] + u[:, Jp1]) / (hθ^2)))

    return rhs
end

function laplace1(u, α, r, hr, hθ)
    rhs = zeros(size(u))
    for i in 2:nr-1
        for j in 2:nθ-1
            rhs[i, j] = α * (((u[i+1, j] - 2 * u[i, j] + u[i-1, j]) / (hr^2)) +
                             (1 / r[i]) * ((u[i+1, j] - u[i-1, j]) / (2 * hr)) +
                             (1 / r[i]^2) * ((u[i, j+1] - 2 * u[i, j] + u[i, j-1]) / (hθ^2)))
        end
    end

    return rhs
end

function update_soln(nt)

    # Numerical solution
    u = copy(u0)
    u_hist = [copy(u)]


    for tstep in 1:nt
        #time = tstep*nt # TODO: implement this in hdf5 file
        # Compute the RHS (2D Laplace)
        rhs = laplace(u, a, r, hr, hθ)

        # Update solution
        u .= u + rhs * ht
        u = apply_BC(u)


        push!(u_hist, copy(u))

        # Stop Criteria
        stop_criteria = sum(rhs[2:end-1, 2:end-1]) # TODO: FINISH THE EXPLANATION: there is a mask since we do not consider the walls of the dommain , see zap!

        if stop_criteria < 1e-3
            println(stop_criteria)
            break
        end

    end

    return u_hist

end
#-----------------------------------------------------------------------------------------------------------------------------------