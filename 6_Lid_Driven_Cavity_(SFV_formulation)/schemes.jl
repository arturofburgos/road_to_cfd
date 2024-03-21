# Here we list the functions needed to solve the numerical problem:
#-----------------------------------------------------------------------------------------------------------------------------------
# Heat Equation 2D:
#=================#
#   Vectorized:   #
#=================#
# For Dirichlet Boundary Conditions 
# and no comparison with the analytical solution
module schemes

function L1norm(new, old)
    norm = sum(abs.(new - old))
    return norm
end



function apply_BC(; psi, w, U, h)

    w[2:end-1, end] .= (-2/h^2)*psi[2:end-1,end-1] .- 2*U/h # psi_top
    w[2:end-1, 1] .= (-2/h^2)*psi[2:end-1,2] # psi_bottom
    w[1, 2:end-1] .= (-2/h^2)*psi[2,2:end-1] # psi_left  
    w[end, 2:end-1] .= (-2/h^2)*psi[end-1,2:end-1] # psi_right

    return w
end


function streamfunction_poisson(psi, w, h)
    psi_old = copy(psi)
    for i in 2:nx-1
        for j in 2:ny-1
            psi[i, j] = 0.25 * (psi[i+1,j] + psi[i-1,j] + psi[i,j+1] + psi[i,j-1] + (h^2)*w[i,j])
        end
    end

    return psi, psi_old
end

function update_soln(psi, w, tn)

    # # Numerical solution
    # u = copy(u0)
    # u_hist = [copy(u)]

    # # # Exact solution
    # # u_exact = copy(u0)
    # # u_hist_exact = [copy(u_exact)]

    for tstep = 1:5
        time = tstep*ht
        psi, psi_old = streamfunction_poisson(psi, w, h)

        w = apply_BC(psi)

        w_old = copy(w)

        for i = 2:n-1
            for  j = 2:n-1
                w[i,j] = alpha*((w[i-1,j] + w[i+1,j] + w[i,j-1] + w[i,j+1]
                                 - 0.25*Re*(psi[i,j+1] - psi[i,j-1])*(w[i+1,j] - w[i-1,j])
                                 + 0.25*Re*(psi[i+1,j] - psi[i-1,j])*(w[i,j+1] - w[i,j-1]))/4)
                                 +(1-alpha)*w[i,j]

            end
        end

        l1_norm_psi = L1norm(psi, psi_old)
        #println("l1_norm_psi = $l1_norm_psi")
        l1_norm_w = L1norm(w, w_old)

        #println("Iteration   = $tn")
        tn += 1


        
    end

    println("Iteration   = $tn")
    println("l1_norm_psi = $l1_norm_psi")
    println("l1_norm_w   = $l1_norm_w")

    return w, psi

end


function update_soln_2(nt)

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

end
#-----------------------------------------------------------------------------------------------------------------------------------