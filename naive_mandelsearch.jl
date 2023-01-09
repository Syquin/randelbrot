include("fractal_img_parallel.jl")
using Images
using ImageIO
using ImageFiltering
using Statistics


"""
mbPointOfInterest(xRes, center, width, its)

Given a window to view and number of iterations in the Mandelbrot set, use the discrete Laplacian to highlight edges. 
Pick the indices of a random point, weighted by the square of the discrete Laplacian.
"""
function mbPointOfInterest(xRes, center, width, its)
    imgmat = graybrotH(xRes, center, width, its)
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

function jlPointOfInterest(xRes, parameter, center, width, its)
    imgmat = grayliaH(xRes, parameter, center, width, its)
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
showmbPointsOfInterest(points, xRes, center, width, its)

Demonstrate mbPointOfInterest by highlighting points number of chosen points in red.
"""
function showmbPointsOfInterest(points, xRes, center, width, its)
    f = t -> t*RGB(1,1,1)
    img = f.(imfilter(graybrotH(xRes, center, width, its), Kernel.Laplacian()))
    for x = 1:points 
        (i, j) = mbPointOfInterest(xRes, center, width, its)
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
showjlPointsOfInterest(points, xRes, center, width, its)

Demonstrate jlPointOfInterest by highlighting points number of chosen points in red.
"""
function showjlPointsOfInterest(points, xRes, parameter, center, width, its)
    f = t -> t*RGB(1,1,1)
    img = f.(imfilter(grayliaH(xRes, parameter, center, width, its), Kernel.Laplacian()))
    for x = 1:points 
        (i, j) = jlPointOfInterest(xRes, parameter, center, width, its)
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
function mbRandZoom(endres, rounds, center, width, zoomfac, its)
    res = zoomfac*10

    center1 = center
    width1 = width 
    for k in 1:rounds 
        dx = width1/res
        topleft = center1 + width1/2*(im-1) 
    
        z(k, l) = topleft + l*dx - k*dx*im
        
        (i, j) = mbPointOfInterest(res, center1, width1, its)
        width1 = width1/zoomfac
        center1 = z(i,j)
    end 
    colorInt = rand(1:5)
    colorList = [sWiki, sRed, sPurple, sGreen, sBlue]
    mandelbrot(endres, center1, width1, its, colorList[colorInt])
end 

function jlRandZoom(endres, rounds, parameter, center, width, zoomfac, its)
    res = zoomfac*10

    center1 = center
    width1 = width 
    for k in 1:rounds 
        dx = width1/res
        topleft = center1 + width1/2*(im-1) 
    
        z(k, l) = topleft + l*dx - k*dx*im
        
        (i, j) = jlPointOfInterest(res, parameter, center1, width1, its)
        width1 = width1/zoomfac
        center1 = z(i,j)
    end 
    colorInt = rand(1:5)
    colorList = [sWiki, sRed, sPurple, sGreen, sBlue]
    juliaSet(endres, parameter, center1, width1, its, colorList[colorInt])
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


#jlRandZoom(500, 4, -1.779 + 0.000im,0im, 0.2, 5, 25000)