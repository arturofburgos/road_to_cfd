function rarefaction(x)
    u = ones(size(x))
    indx1 = findall(x -> x >= 0 && x <= 1, x)
    indx2 = findall(x -> x > 1, x)
    u[indx1] .= 1 .+ x[indx1]
    u[indx2] .= 2
    
    return u
end

function shock(x)
    u = ones(size(x))
    indx1 = findall(x -> x >= 0 && x <= 1, x)
    indx2 = findall(x -> x < 0, x)
    u[indx1] .= 2 .- x[indx1]
    u[indx2] .= 2
    
    return u
end


function initial_condition(x, condition)
    if condition == "shock"
        u = shock(x)

    elseif condition == "rarefaction"
        u = rarefaction(x)

    end
    u0 = copy(u)
    return u, u0
end


function exact_rarefraction(x, time, u0)
    uexact = ones(size(u0))
    for i in eachindex(x)
        if x[i] <= time
            uexact[i] = 1
        elseif (x[i] >= time) && (x[i] <= 1 + 2*time)
            uexact[i] = (1 + x[i]) / (1 + time)
        elseif x[i] >= 1 + 2*time
            uexact[i] = 2
        end
    end
    return uexact
end

function exact_shock(x, time, u0)
    uexact = ones(size(u0))
    for i in eachindex(x)
        if time < 1
            if x[i] < 2*time
                uexact[i] = 2
            elseif (x[i] >= 2*time) && (x[i] <= 1 + time)
                uexact[i] = (2 - x[i]) / (1 - time)
            elseif x[i] > 1 + time
                uexact[i] = 1
            end
        end
        if time > 1
            if x[i] < 3/2*time + 1/2
                uexact[i] = 2
            else
                uexact[i] = 1
            end
        end
    end
    return uexact
end

function f(u)
    return u^2/2
end

function fprime(u)
    return u
end