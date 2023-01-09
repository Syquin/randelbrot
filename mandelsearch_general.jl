include("fractal_img_parallel.jl")
using Images
using ImageIO
using ImageFiltering
using Statistics

function zoomWindow(window::FractalWindow, zoomfactor, targeti, targetj)
    z(i, j) = indexToComplex(window, i, j)
    newCenter = z(targeti,targetj)
    newWidth = window.width/zoomfactor
    FractalWindow(window.xRes, window.yRes, newCenter, newWidth)
end 

"""
imagePointOfInterest(xRes, center, width, its)

Given an image, use the discrete Laplacian to highlight edges. Pick the indices of a random point, weighted by the square of the discrete Laplacian.
"""
function imagePointOfInterest(image)
    (yRes, xRes) = size(image)
    moddedimg = imfilter(image, Kernel.Laplacian())
    total = 0
    for j in 1:xRes 
        for i in 1:yRes 
            total += (moddedimg[i, j])^2
        end
    end # for some reason the built-in summing method over the matrix failed repeatedly :(
    thresh = total*rand()
    s = 0
    for j in 1:xRes 
        for i in 1:yRes 
            s += (moddedimg[i, j])^2
            if s >= thresh
                return (i, j)
            end
        end
    end
end  
"""
showPointsOfInterest(points, image)

Demonstrate imagePointOfInterest by highlighting points number of chosen points in red.
"""
function showPointsOfInterest(points, image)
    f = t -> t*RGB(1,1,1)
    img = imfilter(image, Kernel.Laplacian())
    for x = 1:points 
        (i, j) = imagePointOfInterest(image)
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
mbRandZoom(endres, rounds, center, width, zoomfac, its)

Start with the given center, width, window. Then zoom in with a factor of zoomfac \"rounds\" number of times. Return an image of this location 
with resolution endres x endres.
"""
function randZoom(system, rounds, startWindow, endres, zoomfac, its)
    aspectRatio = startWindow.yRes/startWindow.xRes
    window1 = FractalWindow(endres/10, endres/10, startWindow.center, startWindow.width )
    width1 = startWindow.width 
    center1 = startWindow.center 
    for k in 1:rounds 
        image = prerender(system, window1, its)
        
        (i, j) = imagePointOfInterest(image)
    
        width1 = width1/zoomfac
        center1 = indexToComplex(window1, i, j)
        window1 = FractalWindow(endres/10, endres/10, center1, width1)
    end 
    finalWindow = FractalWindow(endres, endres*aspectRatio, center1, width1)
    colorInt = rand(1:5)
    colorList = [sWiki, sRed, sPurple, sGreen, sBlue]
    render(system, finalWindow, its, colorList[colorInt])
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
function benchmark(n, endres, rounds, center, width, zoomfac, its)
    results = [ @elapsed mbRandZoom(endres, rounds, center, width, zoomfac, its) for i =1:n ]
    println("Average runtime was ", mean(results), "s.")
    println("Standard deviation of runtime was ", std(results), "s.")
end 