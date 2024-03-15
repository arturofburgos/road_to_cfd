# Here we list the functions needed to solve the numerical problem:
#-----------------------------------------------------------------------------------------------------------------------------------
# Godonov:
#=================#
# Non-vectorized: #
#=================#
# For Dirichlet Boundary Conditions 
# and no comparison with the analytical solution

function godonov_flux(u, f)
    uplus = u[Jp1]
    uminus = u[J]

    flux = similar(u)

    for i in 1:length(u)
        if uplus[i] < uminus[i] #rarefaction
            flux[i] = min(f(uminus[i]), f(uplus[i]))
        else
            flux[i] = max(f(uminus[i]), f(uplus[i]))
        end
    end
    
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
        flux = godonov_flux(u,f)

        u[mask] .= (u - ((ht / hx) * (flux[J] - flux[Jm1])))[mask]

        # if condition == "shock"
        #     uexact[mask] = exact_shock(x, time, u0)[mask]

        # elseif condition == "rarefaction"
        #     uexact[mask] = exact_rarefraction(x, time, u0)[mask]
        # end

        push!(u_hist, copy(u))
        # push!(u_hist_exact, copy(u_exact))

    end

    return u_hist, u_hist_exact

end


        






#-----------------------------------------------------------------------------------------------------------------------------------