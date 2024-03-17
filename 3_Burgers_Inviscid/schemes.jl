# Here we list the functions needed to solve the numerical problem:
#-----------------------------------------------------------------------------------------------------------------------------------
# LLF:
#=================#
# Non-vectorized: #
#=================#
# For Dirichlet Boundary Conditions 
# and no comparison with the analytical solution


function LLF(u, f, fprime)

    uplus = u[Jp1]
    uminus = u[J]
    α = max.(abs.(fprime(uplus)),abs.(fprime(uminus)))
    flux = (f.(uplus)+f.(uminus))/2 - α/2 .*(uplus-uminus)

    return flux

end

function update_soln(nt, condition)
    
    # Numerical solution
    u = copy(u0)
    u_hist = [copy(u)]

    # Exact solution
    u_exact = copy(u0)
    u_hist_exact = [copy(u_exact)]

    for tstep in 1:nt
        time = tstep*ht
        flux = LLF(u,f,fprime)

        u[mask] .= (u - ((ht / hx) * (flux[J] - flux[Jm1])))[mask]

        if condition == "shock"
            u_exact[mask] = exact_shock(x, time, u0)[mask]

        elseif condition == "rarefaction"
            u_exact[mask] = exact_rarefraction(x, time, u0)[mask]
        end

        push!(u_hist, copy(u))
        push!(u_hist_exact, copy(u_exact))

    end

    return u_hist, u_hist_exact

end
#-----------------------------------------------------------------------------------------------------------------------------------