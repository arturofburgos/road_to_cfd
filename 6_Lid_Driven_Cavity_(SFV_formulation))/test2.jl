using Plots, LaTeXStrings, Measures # Possible to use PyPlot too.
include("schemes.jl")
include("plotdefaults.jl")


function L1norm(new, old)
    norm = sum(abs.(new - old))
    return norm
end



function calculate_psi_internal(psi, w, h)
    psi_new = copy(psi)
    for i in 2:size(psi, 1)-1
        for j in 2:size(psi, 2)-1
            psi_new[i, j] = 0.25 * ((psi[i + 1, j] + psi[i - 1, j] + psi[i, j + 1] + psi[i, j - 1]) + (h^2) * w[i, j])
        end
    end
    return psi_new
end

function calculate_vorticity_wall(psi, w, h, U)
    w_new = copy(w)
    w_new[2:end-1, 1] = (-2. / (h^2)) * psi[2:end-1, 2]  # bottom
    w_new[2:end-1, end] = (-2. / (h^2)) * psi[2:end-1, end-1] .- 2. * U / h  # top
    w_new[1, 2:end-1] = (-2. / (h^2)) * psi[2, 2:end-1]  # left
    w_new[end, 2:end-1] = (-2. / (h^2)) * psi[end-1, 2:end-1]  # right
    return w_new
end

function calculate_vorticity_internal(psi, w, alpha, Re)
    w_new = copy(w)
    for i in 2:size(w, 1)-1
        for j in 2:size(w, 2)-1
            w_new[i, j] = alpha * ((w[i - 1, j] + w[i + 1, j] + w[i, j - 1] + w[i, j + 1]
                     - 0.25 * Re * (psi[i, j + 1] - psi[i, j - 1]) * (w[i + 1, j] - w[i - 1, j])
                     + 0.25 * Re * (psi[i + 1, j] - psi[i - 1, j]) * (w[i, j + 1] - w[i, j - 1])) / 4) + (1 - alpha) * w[i, j]
        end
    end
    return w_new
end

function main_loop(psi, w, h, U, alpha, Re, l1_target)
    tn = 0
    while true
        psiold = copy(psi)
        psi = calculate_psi_internal(psi, w, h)
        
        w = calculate_vorticity_wall(psi, w, h, U)
        
        wold = copy(w)
        w = calculate_vorticity_internal(psi, w, alpha, Re)
        
        l1_norm_psi = L1norm(psi, psiold)
        l1_norm_omega = L1norm(w, wold)
        tn += 1
        
        if l1_norm_omega ≤ l1_target && l1_norm_psi ≤ l1_target
            break
        end
    end
    return psi, w, tn, l1_norm_psi, l1_norm_omega
end


#-----------------------------------------------------------------------------------------------------------------------------------

#=======================#
# Initialize variables: #
#=======================#
L = 1



N = 5

# Spatial grid for the numerical solution 
x = range(0, L, N)
y = range(0, L, N)
# hx = x[2] - x[1]
# hy = y[2] - y[1]
h = x[2] - x[1]





U = 1
Re = 10
alpha = 1.5
tn = 0 
psi = zeros(N,N)
w = zeros(N,N)

psi0 = copy(psi)
w0 = copy(w)


# Convergence criteria
l1_norm_omega = 1.0
l1_norm_psi = 1.0
l1_target = 1e-6




println("Iteration = ", tn)
println("l1_norm_psi = ", l1_norm_psi)
println("l1_norm_omega = ", l1_norm_omega)

psi, w, tn, l1_norm_psi, l1_norm_omega = main_loop(psi, w, h, U, alpha, Re, l1_target)


# # Call update_solution
# w, psi = update_soln(psi, w, tn)


