include("fractal_img_parallel.jl")
using Images
using ImageIO
using ImageFiltering
using Statistics

"""
graybrotH(xRes, xMin, xMax, yMax, its)

Helper function which generates only the raw numerical data for the mandelbrot set before turning it into an RGB image. 
Used here for edge detection analysis when finding points of interest.
"""
function graybrotH(xRes, xMin, xMax, yMax, its)
    xDist = xMax-xMin

    yRes = floor(Int, xRes)
    dx = xDist/xRes
    x(t) = xMin + (t*dx)
    y(t) = yMax - t*dx

    matrix2 = zeros(Float16, (yRes, xRes))
    Threads.@threads for i in 1:xRes
        Threads.@threads for j in 1:yRes
            matrix2[j, i] = (mandelLeaveTime(x(i), y(j), its, 2))/its
        end
    end
    matrix2
end

"""
smartMandelbrot(xRes, xMin, xMax, yMax, its, coloring = sWiki)

An improvement of my mandelbrot generation function. Previous versions had bad performance when they contained many points of the 
Mandelbrot set itself, due to their maximal number of iterations computed. Here they 
"""
function smartMandelbrot(xRes, xMin, xMax, yMax, its, coloring = sWiki)
    sketch = graybrotH(div(xRes,10), xMin, xMax, yMax, its)
    edginess = abs.(imfilter(sketch, Kernel.Laplacian()))
    marked = zeros(Bool, div(xRes,10),div(xRes,10))
    for i in 1:div(xRes,10)
        for j in 1:div(xRes,10)
            if sketch[i,j] == 0 && edginess[i,j] < 0.01
                marked[i,j] = True 
            end 
        end 
    end 
    xDist = xMax-xMin

    yRes = floor(Int, xRes)
    dx = xDist/xRes

    mark(i,j) = marked[div(i+9, 10),div(j+9, 10)]

    x(t) = xMin + (t*dx)
    y(t) = yMax - t*dx

    matrix2 = zeros(Float16, (yRes, xRes))
    Threads.@threads for i in 1:xRes
        Threads.@threads for j in 1:yRes
            if mark(j, i)
                matrix2[j,i] = 0
            else
                matrix2[j, i] = (mandelLeaveTime(x(i), y(j), its, 2))/its
            end 
        end
    end
    coloring.(matrix2)
end

"""
mbPointOfInterest(xRes, xMin, xMax, yMax, its)

Given a window to view and number of iterations in the Mandelbrot set, use the discrete Laplacian to highlight edges. 
Pick the indices of a random point, weighted by the square of the discrete Laplacian.
"""
function mbPointOfInterest(xRes, xMin, xMax, yMax, its)
    imgmat = graybrotH(xRes, xMin, xMax, yMax, its)
    moddedimg = imfilter(imgmat, Kernel.Laplacian())
    total = 0
    for j in 1:xRes
        for i in 1:xRes
            total += (moddedimg[i, j])^2
        end
    end
    thresh = total*rand()
    s = 0
    for j in 1:xRes
        for i in 1:xRes
            s += (moddedimg[i, j])^2
            if s >= thresh
                return (i, j)
            end
        end
    end
end  
"""
showPointsOfInterest(points, xRes, xMin, xMax, yMax, its)

Demonstrate mbPointOfInterest by highlighting points number of chosen points in red.
"""
function showPointsOfInterest(points, xRes, xMin, xMax, yMax, its)
    f = t -> t*RGB(1,1,1)
    img = f.(imfilter(graybrotH(xRes, xMin, xMax, yMax, its), Kernel.Laplacian()))
    for x = 1:points 
        (i, j) = mbPointOfInterest(xRes, xMin, xMax, yMax, its)
        for k = (i-2):(i+2)
            for l = (j-2):(j+2)
                if (0 < k) && (k<=xRes) &&(0 < l) && (l<=xRes)
                    img[k, l] = RGB(1, 0, 0)
                end
            end
        end
    end
    img
end 

"""
mbRandZoom(endres, rounds, xMin, xMax, yMax, zoomfac, its)

Start with the given xMin, xMax, yMax window. Then zoom in with a factor of zoomfac \"rounds\" number of times. Return an image of this location 
with resolution endres x endres.
"""
function mbRandZoom(endres, rounds, xMin, xMax, yMax, zoomfac, its)
    res = zoomfac*10

    xMin1 = xMin
    xMax1 = xMax
    yMax1 = yMax
    for k in 1:rounds 
        xDist = xMax1-xMin1
        dx = xDist/res
        x(t) = xMin1 + (t*dx)
        y(t) = yMax1 - t*dx
        (i, j) = mbPointOfInterest(res, xMin1, xMax1, yMax1, its)
        (a, b) = (x(j), y(i))
        xDist = xDist/zoomfac
        xMin1 = a - xDist/2 
        xMax1 = a + xDist/2
        yMax1 = b + xDist/2
    end 
    colorInt = rand(1:5)
    colorList = [sWiki, sRed, sPurple, sGreen, sBlue]
    smartMandelbrot(endres, xMin1, xMax1, yMax1, its, colorList[colorInt])
end 

"""
quick helper function to save 1-1000 files automatically
"""
function threeDigiter(n)
    if n < 10
        return "00"*string(n)
    elseif n<100
        return "0"*string(n)
    else
        return string(n)
    end
end
"""
quick test for speed
"""
function benchmark(n, endres, rounds, xMin, xMax, yMax, zoomfac, its)
    results = [ @elapsed mbRandZoom(endres, rounds, xMin, xMax, yMax, zoomfac, its) for i =1:n ]
    println("Average runtime was ", mean(results), "s.")
    println("Standard deviation of runtime was ", std(results), "s.")
end 

benchmark(100, 500, 4, -1.4,-1.3, 0.1, 10, 25000)