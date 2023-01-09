include("sGradients.jl")
using ImageFiltering

"""
mandelbrot defining function
"""
function mbf(z::Complex, c::Complex) 
    z^2 + c
end 

"""
burning ship defining function
"""
function bsf(z::Complex, c::Complex) 
    (x,y) =(real(z), imag(z))
    (abs(x) + abs(y)*im)^2 + c 
end 


 #=       Fix a point c on the complex plane. Let f_c be the function 
    z -> z^2 + c. The number c is in the Mandelbrot set if (f_c)^n(0) 
    does not go to infinity. The below function checks if the orbit 
    goes to infinity by looking up to n steps. If it diverges, it returns
    the number of steps it takes for the point 0 to have an absolute value above 2. 
    The reason to return this value is to color the graph later.
=#


function mandelLeaveTime(f, c, n, bound)
    z = 0 + 0*im
    b2 = bound^2
    if abs2(z) >= b2
        return 0
    end
    for i in 1:n
        z = f(z,c)
        if abs2(z) >=  b2
            return i
        end
    end
    return n+1
end

# some functions which make color gradients work
#makeColorTracker(n) = t -> ((n+ 1 -t)/(n+1))^
makeColorTracker(n) = t ->  t < 0 ? 0 : (t == n+1 ? 0 : ((t^(1/5)) % 10)/10)


#=
        The Julia set at a point c is defined similarly. Define f_c as before.
    The Julia set consists of the points z such that (f_c)^n(z) does not go
    to infinity. The value returned is decided in the same way as before.
=#
function juliaLeaveTime(f, z, c, n, bound)
    b2 = bound^2
    if abs2(z) >=  b2
        return 0
    end
    z1 = z
    for i in 1:n
        z1 = f(z1,c)
        if abs(z1) >=  b2
            return i
        end
    end
    return n+1
end


"""
graybrotH(xRes, xMin, xMax, yMax, its)

Helper function which generates only the raw numerical data for the mandelbrot set before turning it into an RGB image. 
Used here for edge detection analysis when finding points of interest.
"""
function graybrotH(xRes, center, width, its)
    dx = width/xRes
    topleft = center + width/2*(im-1) 

    z(k, l) = topleft + l*dx - k*dx*im

    matrix2 = zeros(Float16, (xRes, xRes))
    Threads.@threads for i in 1:xRes
        Threads.@threads for j in 1:xRes
            matrix2[i, j] = (mandelLeaveTime(mbf, z(i,j), its, 2))/its
        end
    end
    matrix2
end

"""
grayliaH(xRes, xMin, xMax, yMax, its)

Helper function which generates only the raw numerical data for the mandelbrot set before turning it into an RGB image. 
Used here for edge detection analysis when finding points of interest.
"""
function grayliaH(xRes, parameter, center, width, its)
    dx = width/xRes
    topleft = center + width/2*(im-1) 

    z(k, l) = topleft + l*dx - k*dx*im

    matrix2 = zeros(Float16, (xRes, xRes))
    Threads.@threads for i in 1:xRes
        Threads.@threads for j in 1:xRes
            matrix2[i, j] = (juliaLeaveTime(mbf, z(i,j), parameter, its, 2))/its
        end
    end
    matrix2
end



function mandelbrotH(f, xRes, center, width, its)
    dx = width/xRes
    topleft = center + width/2*(im-1) 
    

    sketch = graybrotH(div(xRes,10), center, width, its)
    edginess = abs.(imfilter(sketch, Kernel.Laplacian()))
    marked = zeros(Bool, div(xRes,10),div(xRes,10))

    for i in 1:div(xRes,10)
        for j in 1:div(xRes,10)
            if sketch[i,j] == 0 && edginess[i,j] < 0.01
                marked[i,j] = True 
            end 
        end 
    end 
    
    z(k, l) = topleft + l*dx - k*dx*im


    mark(i,j) = marked[div(i+9, 10),div(j+9, 10)]
    c = makeColorTracker(its)
    matrix2 = zeros(Float16, (xRes, xRes))
    Threads.@threads for i in 1:xRes
        Threads.@threads for j in 1:xRes
            if mark(i, j)
                matrix2[i,j] = 0
            else
                matrix2[i, j] = c(mandelLeaveTime(f, z(i,j), its, 2))
            end 
        end
    end
    matrix2
end
mandelbrot(xRes, center, width, its, coloring=sThermal) = coloring.(mandelbrotH(mbf, xRes, center, width, its))


# Generate a graph of the Burning Ship with x values in the range xr,
# y values in the range yr, and a maximum number of n iterations.

burningship(xRes, center, width, its, coloring=sThermal) = coloring.(mandelbrotH(bsf, xRes, center, width, its))


# Generate a graph of the Julia set with c = (real + i*imaginary),
# x values in the range xr, y values in the range yr, and a maximum number of n iterations.

function juliaSetH(xRes, c, center, width, its)
    dx = width/xRes
    topleft = center + width/2*(im-1) 
    

    sketch = graybrotH(div(xRes,10), center, width, its)
    edginess = abs.(imfilter(sketch, Kernel.Laplacian()))
    marked = zeros(Bool, div(xRes,10),div(xRes,10))

    for i in 1:div(xRes,10)
        for j in 1:div(xRes,10)
            if sketch[i,j] == 0 && edginess[i,j] < 0.01
                marked[i,j] = True 
            end 
        end 
    end 
    
    z(k, l) = topleft + l*dx - k*dx*im


    mark(i,j) = marked[div(i+9, 10),div(j+9, 10)]
    col = makeColorTracker(its)
    matrix2 = zeros(Float16, (xRes, xRes))
    Threads.@threads for i in 1:xRes
        Threads.@threads for j in 1:xRes
            if mark(i, j)
                matrix2[i,j] = 0
            else
                matrix2[i, j] = col(juliaLeaveTime(mbf, z(i,j), c, its, 2))
            end 
        end
    end
    matrix2
end
juliaSet(xRes, c, center, width, its, coloring=sWiki) = coloring.(juliaSetH(xRes, c, center, width, its))