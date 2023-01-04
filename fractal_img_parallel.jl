include("sGradients.jl")
              

 #=       Fix a point c on the complex plane. Let f_c be the function 
    z -> z^2 + c. The number c is in the Mandelbrot set if (f_c)^n(0) 
    does not go to infinity. The below function checks if the orbit 
    goes to infinity by looking up to n steps. If it diverges, it returns
    the number of steps it takes for the point 0 to have an absolute value above 2. 
    The reason to return this value is to color the graph later.
=#

function mandelLeaveTime(x0, y0, n, bound)
    x1 = 0
    y1 = 0
    b2 = bound^2
    if x0^2 + y0^2 >=  b2
        return 0
    end
    for i in 1:n
        xtem = (x1 + y1)*(x1-y1) + x0
        ytem = 2*x1*y1 + y0
        x1 = xtem
        y1 = ytem
        if x1^2 + y1^2 >=  b2
            return i
        end
    end
    return n+1
end

# some functions which make color gradients work
#makeColorTracker(n) = t -> ((n+ 1 -t)/(n+1))^
makeColorTracker(n) = t ->  t < 0 ? 0 : (t == n+1 ? 0 : (sqrt(t) % 17)/17)
makeColorTracker2(cyc, n) = t -> ( t == n+1 ? 0 : ((t%cyc)/cyc))

makeColorTrackerJ(n) = t ->  t < 0 ? 0 : (t == n+1 ? 0 : sqrt(t) % 7)/7


#=
        The Julia set at a point c is defined similarly. Define f_c as before.
    The Julia set consists of the points z such that (f_c)^n(z) does not go
    to infinity. The value returned is decided in the same way as before.
=#
function juliaLeaveTime(x, y, x0, y0, n, bound)
    x1 = x
    y1 = y
    b2 = bound^2
    if x0^2 + y0^2 >=  b2
        return 0
    end
    for i in 1:n
        xtem = (x1 + y1)*(x1-y1) + x0
        ytem = 2*x1*y1 + y0
        x1 = xtem
        y1 = ytem
        if x1^2 + y1^2 >=  b2
            return i
        end
    end
    return n+1
end



#=
    The Burning Ship fractal is defined similarly to the Mandelbrot set, except that we 
take Mandel(x0, y0, |x|, -|y|) instead. See this link for more info.

https://en.wikipedia.org/wiki/Burning_Ship_fractal

It's my personal favorite fractal as of right now.
=#
bShip(x0, y0, x, y) = Mandel(x0, y0, abs(x), -abs(y)) 
function bShipLeaveTime(x0, y0, n, bound)
    x1 = 0
    y1 = 0
    b2 = bound^2
    if x0^2 + y0^2 >=  b2
        return 0
    end
    for i in 1:n
        xtem = (x1 + y1)*(x1-y1) + x0
        ytem = 2*x1*y1 - y0
        x1 = abs(xtem)
        y1 = abs(ytem)
        if x1^2 + y1^2 >=  b2
            return i
        end
    end
    return n+1
end



# Generate a graph of the mandelbrot set with x values in the range xr,
# y values in the range yr, and a maximum number of n iterations.
using LoopVectorization


function mandelbrotH(xRes, xMin, xMax, yMax, its)
    xDist = xMax-xMin
    c = makeColorTracker(its)

    yRes = floor(Int, xRes/4*3)
    dx = xDist/xRes
    x(t) = xMin + (t*dx)
    y(t) = yMax - t*dx

    matrix2 = zeros(Float64, (yRes, xRes))
    Threads.@threads for i in 1:xRes
        Threads.@threads for j in 1:yRes
            matrix2[j, i] = c(mandelLeaveTime(x(i), y(j), its, 2))
        end
    end
    matrix2
end
mandelbrot(xRes, xMin, xMax, yMax, its, coloring=sThermal) = coloring.(mandelbrotH(xRes, xMin, xMax, yMax, its))


# Generate a graph of the Burning Ship with x values in the range xr,
# y values in the range yr, and a maximum number of n iterations.

function burningShipH(xRes, xMin, xMax, yMax, its)
    xDist = xMax-xMin
    c = makeColorTracker(its)

    yRes = floor(Int, xRes/4*3)
    dx = xDist/xRes
    x(t) = xMin + (t*dx)
    y(t) = yMax - t*dx

    matrix2 = zeros(Float64, (yRes, xRes))
    Threads.@threads for i in 1:xRes
        Threads.@threads for j in 1:yRes
            matrix2[j, i] = c(bShipLeaveTime(x(i), y(j), its, 2))
        end
    end
    matrix2
end
burningShip(xRes, xMin, xMax, yMax, its, coloring=sThermal) = coloring.(burningShipH(xRes, xMin, xMax, yMax, its))


# Generate a graph of the Julia set with c = (real + i*imaginary),
# x values in the range xr, y values in the range yr, and a maximum number of n iterations.

function juliaSetH(real, imaginary, xRes, xMin, xMax, yMax, its)
    xDist = xMax-xMin
    c = makeColorTrackerJ(its)

    yRes = floor(Int, xRes/4*3)
    dx = xDist/xRes
    x(t) = xMin + (t*dx)
    y(t) = yMax - t*dx

    matrix2 = zeros(Float64, (yRes, xRes))
    Threads.@threads for i in 1:xRes
        Threads.@threads for j in 1:yRes
            matrix2[j, i] = c(juliaLeaveTime(x(i), y(j), real, imaginary,  its, 2))
        end
    end
    matrix2
end
juliaSet(real, imaginary, xRes, xMin, xMax, yMax, its, coloring=sThermal) = coloring.(juliaSetH(real, imaginary, xRes, xMin, xMax, yMax, its))
