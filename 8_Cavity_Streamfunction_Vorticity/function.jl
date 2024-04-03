module soln
using Printf

function L1norm(new, old)
    norm = sum(abs.(new - old))
    return norm
end

function apply_BC(; w, psi, U, h)

    w[2:end-1, end] .= (-2 / h^2) * psi[2:end-1, end-1] .- 2 * U / h # psi_top
    w[2:end-1, 1] .= (-2 / h^2) * psi[2:end-1, 2] # psi_bottom
    w[1, 2:end-1] .= (-2 / h^2) * psi[2, 2:end-1] # psi_left  
    w[end, 2:end-1] .= (-2 / h^2) * psi[end-1, 2:end-1] # psi_right

    return w
end

function streamfunction_poisson(; psi, w, h, nx, ny)
    psi_old = copy(psi)
    for i in 2:nx-1
        for j in 2:ny-1
            psi[i, j] = 0.25 * (psi[i+1, j] + psi[i-1, j] + psi[i, j+1] + psi[i, j-1] + (h^2) * w[i, j])
        end
    end

    return psi, psi_old
end

function update_soln(; nx, ny, psi, w, h, N, U, Re, alpha, w_hist, psi_hist, dt ,v)

    l1_norm_omega = Inf
    l1_norm_psi = Inf
    l1_target = 1e-6
    tn = 0

    push!(w_hist, copy(w))
    push!(psi_hist, copy(psi))

    for step = 1:500 # Uncoment if you want to use for while loop, make sure to comment the below one

        # while (l1_norm_omega > l1_target) || (l1_norm_psi > l1_target)
        # Calculate psi at internal grid points (BCs for psi are automatically enforced)
        psi, psiold = streamfunction_poisson(; psi, w, h, nx, ny)


        # psiold = copy(psi)
        # psi = calculate_psi_internal(psi0, w0, h)


        # Calculate Vorticity at the wall; omit corner points due to singularity
        w = apply_BC(; w, psi, U, h)
        wold = copy(w)

        # Calculate Vorticity at internal grid points        
        for i in 2:N-1
            for j in 2:N-1
                w[i, j] = w[i,j] 
                - 
                dt*((psi[i, j+1] - psi[i, j-1])/(2 * h))*((w[i+1, j] -  w[i-1, j])/2*h) 
                - 
                dt*((psi[i+1, j] - psi[i-1, j])/(2 * h))*((w[i, j+1] -  w[i, j-1])/2*h)
                +
                dt*v*((w[i-1, j] + w[i+1, j] + w[i, j-1] + w[i, j+1] - 4*w[i, j])/(h^2))
            end
        end


        # for i in 2:N-1
        #     for j in 2:N-1
        #         w[i, j] = alpha * ((w[i-1, j] + w[i+1, j] + w[i, j-1] + w[i, j+1]
        #                             -
        #                             0.25 * Re * (psi[i, j+1] - psi[i, j-1]) * (w[i+1, j] - w[i-1, j])
        #                             +
        #                             0.25 * Re * (psi[i+1, j] - psi[i-1, j]) * (w[i, j+1] - w[i, j-1])) / 4) + (1 - alpha) * w[i, j]
        #     end
        # end


        l1_norm_psi = L1norm(psi, psiold)
        l1_norm_omega = L1norm(w, wold)
        tn = step + 1
        # tn += 1 # Uncoment if you want to use for while loop, make sure to comment the above one
        push!(w_hist, copy(w))
        push!(psi_hist, copy(psi))

        @printf(
            stderr,
            "iter %6d | l1_norm_psi %.3e | l1_norm_omega %.3e\n",
            tn, l1_norm_psi, l1_norm_omega
        )
    end
    return psi_hist, w_hist, tn
end

export update_soln


end