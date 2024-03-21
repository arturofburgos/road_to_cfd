using Plots, LaTeXStrings, Measures # Possible to use PyPlot too.
include("schemes.jl")
include("plotdefaults.jl")






function apply_BC(psi)

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



function update_soln1(psi, w, tn)

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


function calculate_psi_internal(psi, w, h)
    psi_new = copy(psi)
    for i in 2:size(psi, 1)-1
        for j in 2:size(psi, 2)-1
            psi_new[i, j] = 0.25 * ((psi[i + 1, j] + psi[i - 1, j] + psi[i, j + 1] + psi[i, j - 1]) + (h^2) * w[i, j])
        end
    end
    return psi_new
end


#-----------------------------------------------------------------------------------------------------------------------------------

#=======================#
# Initialize variables: #
#=======================#




L = 1

N = 31


x = range(0, L, N)
y = range(0, L, N)

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

w_hist = []
psi_hist = []


module soln

function L1norm(new, old)
    norm = sum(abs.(new - old))
    return norm
end

function update_soln(;psi, w, h, N, U, Re, alpha, w_hist, psi_hist)
    
    l1_norm_omega = Inf
    l1_norm_psi = Inf
    l1_target = 1e-6
    tn = 0 
    
    push!(w_hist, copy(w))
    push!(psi_hist, copy(psi))

    for step = 1:500 # Uncoment if you want to use for while loop, make sure to comment the below one
        
    # while l1_norm_omega > l1_target || l1_norm_psi > l1_target
        # Calculate psi at internal grid points (BCs for psi are automatically enforced)
        psiold = copy(psi)
        for i in 2:N-1
            for j in 2:N-1
                psi[i, j] = 0.25 * ((psi[i + 1, j] + psi[i - 1, j] + psi[i, j + 1] + psi[i, j - 1]) + (h^2) * w[i, j])
            end
        end


        # psiold = copy(psi)
        # psi = calculate_psi_internal(psi0, w0, h)


        # Calculate Vorticity at the wall; omit corner points due to singularity
        w[2:end-1, 1] = (-2. / (h^2)) * psi[2:end-1, 2]  # bottom
        w[2:end-1, end] = (-2. / (h^2)) * psi[2:end-1, end-1] .- 2. * U / h  # top
        w[1, 2:end-1] = (-2. / (h^2)) * psi[2, 2:end-1]  # left
        w[end, 2:end-1] = (-2. / (h^2)) * psi[end-1, 2:end-1]  # right

        # Calculate Vorticity at internal grid points
        wold = copy(w)
        for i in 2:N-1
            for j in 2:N-1
                w[i, j] = alpha * ((w[i - 1, j] + w[i + 1, j] + w[i, j - 1] + w[i, j + 1]
                        - 0.25 * Re * (psi[i, j + 1] - psi[i, j - 1]) * (w[i + 1, j] - w[i - 1, j])
                        + 0.25 * Re * (psi[i + 1, j] - psi[i - 1, j]) * (w[i, j + 1] - w[i, j - 1])) / 4) + (1 - alpha) * w[i, j]
            end
        end

        l1_norm_psi = L1norm(psi, psiold)
        l1_norm_omega = L1norm(w, wold)
        tn = step +1
        # tn += 1 # Uncoment if you want to use for while loop, make sure to comment the above one
        push!(w_hist, copy(w))
        push!(psi_hist, copy(psi))
    end
    return  psi_hist, w_hist, tn
end
end

println("Iteration = ", tn)
println("l1_norm_psi = ", l1_norm_psi)
println("l1_norm_omega = ", l1_norm_omega)

begin
    psi_hist, w_hist, tn  = soln.update_soln(;psi, w, h, N, U, Re, alpha, w_hist, psi_hist)

end

# Update plot functions --> to be used in the .gif
function update_plot_physical(n, u_hist, x, y)
    contourf(x, y, transpose(u_hist[n]),
        color=:jet, levels=20, aspect_ratio=:equal)
    xlabel!("x")
    ylabel!("y")
end



# Create the .gif's
anim1 = @animate for n in 1:10:length(w_hist)
    update_plot_physical(n, w_hist, x, y)
end

anim2 = @animate for n in 1:10:length(w_hist)
    update_plot_physical(n, psi_hist, x, y)
end


# Save the .gif
gif(anim1, fps=10)

gif(anim2, fps=10)
