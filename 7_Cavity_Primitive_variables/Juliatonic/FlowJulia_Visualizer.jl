using CairoMakie, FileIO, Printf

#### Simulation inputs
length = 4.0
breadth = 4.0
colpts = 111
rowpts = 111


using Printf
using DelimitedFiles

# Go to the Result directory
cwdir = pwd()
dir_path = joinpath(cwdir, "7_Cavity_Primitive_variables/Juliatonic/Result/Snapshots")
# cd(dir_path)

# Go through files in the directory and store filenames
filenames = []
iterations = []

for (root, dirs, files) in walkdir(dir_path)
    for datafile in files
        if occursin("PUV", datafile)
            push!(filenames, datafile)
            no_ext_file = replace(datafile, ".txt" => "") |> strip
            iter_no = parse(Int, split(no_ext_file, "V")[end])
            push!(iterations, iter_no)
        end
    end
end

# Discern the final iteration and interval
initial_iter = minimum(iterations)
final_iter = maximum(iterations)
inter = (final_iter - initial_iter) / (size(iterations)[1] - 1)
number_of_frames = size(iterations)[1]
sorted_iterations = sort(iterations)

# Output the results
@printf("Initial Iteration: %d\n", initial_iter)
@printf("Final Iteration: %d\n", final_iter)
@printf("Interval: %.2f\n", inter)
@printf("Number of Frames: %d\n", number_of_frames)
@printf("Sorted Iterations: %s\n", join(sorted_iterations, ", "))

function ReadDataFile(iteration::Int)
    # Set filename and path according to given iteration
    filename = "PUV$iteration.txt"
    filepath = joinpath(dir_path, filename)
    
    # Load text file as array
    arr = readdlm(filepath, '\t')
    rows, cols = size(arr)
    
    # Define empty arrays for pressure and velocities
    p_p = zeros(rowpts, colpts)
    u_p = zeros(rowpts, colpts)
    v_p = zeros(rowpts, colpts)
    
    # Organize imported array into variables
    p_arr = arr[:, 1]
    u_arr = arr[:, 2]
    v_arr = arr[:, 3]
    
    # Reshape 1D data into 2D
    p_p = reshape(p_arr, (rowpts, colpts))
    u_p = reshape(u_arr, (rowpts, colpts))
    v_p = reshape(v_arr, (rowpts, colpts))
    
    return p_p, u_p, v_p
end


x = range(0, length, colpts)
y = range(0, breadth, colpts)

# p_p, u_p, v_p = ReadDataFile(4000)




function VectorField(;length, breadth, u_p, v_p)
    # Scale p values to range from 0 to lenght_value and from 0 to breadth_value
    function (p) 
    
        scaled_p1 = p[1] * colpts
        scaled_p2 = p[2] * rowpts
        # Clamp scaled values to ensure they stay within valid range
        i = clamp(round(Int, scaled_p1/length), 1, colpts)
        j = clamp(round(Int, scaled_p2/breadth), 1, rowpts)
        
        # Fetch the vector components from U and V matrices
        return Point2f(u_p[i, j], v_p[i, j])

    end
end





# we use `record` to show the resulting video in the docs.
# If one doesn't need to record a video, a normal loop works as well.
# Just don't forget to call `display(fig)` before the loop
# and without record, one needs to insert a yield to yield to the render task

xinterval_vec_field = 0..length
yinterval_vec_field = 0..breadth

x = range(0, length, colpts)
y = range(0, breadth, colpts)

fig = Figure()
ax = Axis(fig[1, 1])

i_start = 2
it = sorted_iterations[i_start]

p_p, u_p, v_p = ReadDataFile(it)
f = Observable(VectorField(; length, breadth, u_p, v_p))

cont = contourf!(ax, x, y, p_p, levels=10, colorrange = (0.001,0.6))
sp = streamplot!(ax, f, xinterval_vec_field, yinterval_vec_field, colormap=:magma)

Colorbar(fig[1, 2], cont)
Colorbar(fig[1, 3], sp)

output_video = "results.gif"
record(fig, output_video, i_start:size(sorted_iterations)[1]; framerate = 10) do i
    it = sorted_iterations[i]
    p_p, u_p, v_p = ReadDataFile(it)
    f[] = VectorField(; length, breadth, u_p, v_p)
    cont[3] = p_p
end



# using CairoMakie

# v(x::Point2{T}, t) where T = Point2{T}(one(T) * x[2] * t, 4 * x[1])

# vt(it) = let 
    
#     p_p, u_p, v_p = ReadDataFile(it)
#     VectorField(; length, breadth, u_p, v_p)
    
# end

# sf = Observable{Any}(vt(sorted_iterations[1]))
# title_str = Observable("t = 0.00")
# sp = streamplot(sf, -2..2, -2..2;
#     linewidth=2, colormap=:magma, axis=(;title=title_str))

# record(sp, "output.mp4", 1:size(sorted_iterations)[1]; framerate=10) do t
#     sf[] = vt(sorted_iterations[t])
#     title_str[] = "t = $(round(t; sigdigits=2))"
# end