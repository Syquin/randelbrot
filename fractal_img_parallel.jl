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
    (abs(x) - abs(y)*im)^2 + c 
end 


 #=       Fix a point c on the complex plane. Let f_c be the function 
    z -> z^2 + c. The number c is in the Mandelbrot set if (f_c)^n(0) 
    does not go to infinity. The below function checks if the orbit 
    goes to infinity by looking up to n steps. If it diverges, it returns
    the number of steps it takes for the point 0 to have an absolute value above 2. 
    The reason to return this value is to color the graph later.
=#


struct FractalWindow  
    xRes :: Int # horizonal resolution of image 
    yRes :: Int # vertical resolution of image 
    center :: Complex # center of image in complex plane 
    width :: Real # width of image in complex plane 
end 


struct FractalSystem 
    f # function CxC to C defining the dynamical system 
    isJulia :: Bool # is it a Julia render rather than Mandelbrot 
    parameter:: Complex # only relevant for Julia set 
end 

FSystem(f, i, c=0im) = FractalSystem(f, i, c)

const Mandelbrot = FSystem(mbf, false)
const BurningShip = FSystem(bsf, false)

function julia(c) 
    FSystem(mbf, true, c)
end 

function burningJulia(c) 
    FSystem(bsf, true, c)
end 

"""
indexToComplex(wind::FractalWindow, i, j)

Given pixel indices i,j and a window wind, get complex position of point.
"""
function indexToComplex(window::FractalWindow, i, j)
    dx = window.width/window.xRes
    aspectRatio = window.yRes/window.xRes 
    height = window.width*aspectRatio
    topleft = window.center - window.width/2 + im*height/2 
    topleft + j*dx - i*dx*im
end 

"""
mandelLeaveTime(f, c, n, bound)

Get the number of iterations to leave the bounded region, trying n iterations, for
a Mandelbrot-like function f, location c.
"""
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

# a map from {1, ..., n} to [0,1] to make color gradients pretty 
# across different scales. something of a placeholder.
makeColorTracker(n) = t ->  t < 0 ? 0 : (t == n+1 ? 0 : (((t^0.25)*log(t+1)) % 7)/7)


#=
        The Julia set at a point c is defined similarly. Define f_c as before.
    The Julia set consists of the points z such that (f_c)^n(z) does not go
    to infinity. The value returned is decided in the same way as before.
=#
"""
juliaLeaveTime(f, z, c, n, bound)

Get the number of iterations to leave the bounded region, trying n iterations, for
a julia set function f, location z, parameter c.
"""
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
graybrotH(window, its)

Helper function which generates only the raw numerical data for the mandelbrot set before turning it into an RGB image. 
Used here for edge detection analysis when finding points of interest.
"""
function graybrotH(window, its) 

    z(k, l) = indexToComplex(window, k, l)

    matrix2 = zeros(Float16, (window.yRes, window.xRes))
    Threads.@threads for i in 1:window.yRes
        Threads.@threads for j in 1:window.xRes
            matrix2[i, j] = (mandelLeaveTime(mbf, z(i,j), its, 2))/its
        end
    end
    matrix2
end

"""
grayliaH(window, its)

Helper function which generates only the raw numerical data for the mandelbrot set before turning it into an RGB image. 
Used here for edge detection analysis when finding points of interest.
"""
function grayliaH(parameter, window, its)
    z(k, l) = indexToComplex(window, k, l)

    matrix2 = zeros(Float16, (window.yRes, window.xRes))
    Threads.@threads for i in 1:window.yRes
        Threads.@threads for j in 1:window.xRes
            matrix2[i, j] = (juliaLeaveTime(mbf, z(i,j), parameter, its, 2))/its
        end
    end
    matrix2
end


"""
mandelbrot helper function. does most of the actual work
"""
function mandelbrotH(f, window, its, bounds = 2)

    sketchWindow = FractalWindow(div(window.xRes,10), div(window.yRes,10), window.center, window.width )
    sketch = graybrotH(sketchWindow, its)

    edginess = abs.(imfilter(sketch, Kernel.Laplacian()))
    marked = zeros(Bool, sketchWindow.yRes, sketchWindow.xRes)

    for i in 1:sketchWindow.yRes
        for j in 1:sketchWindow.xRes
            if sketch[i,j] == 0 && edginess[i,j] < 0.01
                marked[i,j] = true 
            end 
        end 
    end 

    z(k, l) = indexToComplex(window, k, l)

    mark(i,j) = marked[div(i+9, 10),div(j+9, 10)]
    c = makeColorTracker(its)
    matrix2 = zeros(Float16, (window.yRes, window.xRes))
    Threads.@threads for i in 1:window.yRes
        Threads.@threads for j in 1:window.xRes
            if mark(i, j)
                matrix2[i,j] = 0
            else
                matrix2[i, j] = c(mandelLeaveTime(f, z(i,j), its, bounds))
            end 
        end
    end
    matrix2
end
"""
mandelbrot(window, its, coloring=sWiki)

Get a Mandelbrot image with given window and iterations. Can pick color gradient with coloring.
"""
mandelbrot(window, its, coloring=sWiki) = coloring.(mandelbrotH(mbf, window, its))


# Generate a graph of the Burning Ship with x values in the range xr,
# y values in the range yr, and a maximum number of n iterations.

"""
burningShip(window, its, coloring=sWiki)

Get a Burning Ship image with given window and iterations Can pick color gradient with coloring.
"""
burningShip(window, its, coloring=sWiki) = coloring.(mandelbrotH(bsf, window, its))


# Generate a graph of the Julia set with c = (real + i*imaginary),
# x values in the range xr, y values in the range yr, and a maximum number of n iterations.

function juliaSetH(f, c, window, its, bound)

    sketchWindow = FractalWindow(div(window.xRes,10), div(window.yRes,10), window.center, window.width )
    sketch = grayliaH(c, sketchWindow, its)

    edginess = abs.(imfilter(sketch, Kernel.Laplacian()))
    marked = zeros(Bool, sketchWindow.yRes, sketchWindow.xRes)

    for i in 1:sketchWindow.yRes
        for j in 1:sketchWindow.xRes
            if sketch[i,j] == 0 && edginess[i,j] < 0.01
                marked[i,j] = true 
            end 
        end 
    end 
    
    z(k, l) = indexToComplex(window, k, l)

    mark(i,j) = marked[div(i+9, 10),div(j+9, 10)]

    col = makeColorTracker(its)
    matrix2 = zeros(Float16, (window.yRes, window.xRes))
    Threads.@threads for i in 1:window.yRes
        Threads.@threads for j in 1:window.xRes
            if mark(i, j)
                matrix2[i,j] = 0
            else
                matrix2[i, j] = col(juliaLeaveTime(f, z(i,j), c, its, bound))
            end 
        end
    end
    matrix2
end

"""
juliaSet(c, window, its, coloring=sWiki)

Get a julia set image with parameter c, given window, given number of iterations. Can pick color gradient with coloring.
"""
juliaSet(c, window, its, coloring=sWiki) = coloring.(juliaSetH(mbf, c, window, its,2))

function prerender(system, window, its, bounds= 2)
    if system.isJulia
        return juliaSetH(system.f, system.parameter, window, its, bounds)
    else 
        return mandelbrotH(system.f, window, its, bounds)
    end 
end 

render(system, window, its, coloring = sWiki, bounds = 2) = coloring.(prerender(system, window, its, bounds))

