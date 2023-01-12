include("mandelsearch_general.jl")

#check if boolean list's true values are connected 
function isConnected(list)
    firstComponentStarted = false 
    firstComponentEnded = false 
    for b in list
        if b && !firstComponentStarted 
            firstComponentStarted = true 
        elseif !b && firstComponentStarted && !firstComponentEnded
            firstComponentEnded = true 
        elseif b && firstComponentEnded
            return false 
        end 
    end 
    true 
end 

#turn precolored fractal image into boolean based on threshold 
function toThreshold(preRendered, threshold)
    makescut(t) = abs(t) >= threshold 
    makescut.(preRendered)
end 

function demoThresholds(image)
    save("threshExamples\\original.png",map(clamp01nan, sWiki.(image)))
    image = imfilter(image, Kernel.Laplacian())
    save("threshExamples\\originaledge.png",map(clamp01nan, Gray.(abs.(image))))
    color(bool) = bool ? Gray(1) : Gray(0)
    for i = 0:0.05:3
        save("threshExamples\\0."*threeDigiter(floor(Int,i*100))*".png", toThreshold(image, i))
    end 
end 

LRKernel = [-1 0 1;
            -1 0 1;
            -1 0 1]
UDKernel = transpose(LRKernel)

function isMoreHorizontal(preimg)
    image1 = toThreshold(
        imfilter(preimg, Kernel.gaussian(10)),
        0.6
    )
    color(bool) = bool ? 1 : 0
    image = color.(image1)
    horimat = imfilter(image, LRKernel)
    vertmat = imfilter(image, UDKernel)
    leftness = sum(horimat .* horimat)
    upness = sum(vertmat .* vertmat)
    leftness >= upness
end 

function getPointsOfInterest(points, image)
    [imagePointOfInterest(image) for i = 1:points ]
end 

function interpolatePointsOfInterest(points, degree, system, window, its)
    image = graybrotH(system, window, its)
    inds = getPointsOfInterest(points, image)
    z(i, j) = indexToComplex(window, i, j)
    parts(comp, i) = (i == 1) ? real(comp) : imag(comp)
    interpolator = [ inds[i][3-j] for i = 1:points, j = 1:2]
    vec = polyRegress(interpolator, degree)
    f = polyOffset( vec , 0)

    (yRes, xRes) = size(image)
    xCent = real(window.center)
    dx = window.width/xRes
    xInd(x) = (x-xCent)/dx + xRes/2

    
    func(x) = imag(z(f(xInd(x)),xInd(x)))  # not the most optimized but i'm lazy today
end 

function demonstrateInterpolate(points, degree, system, window, its)
    # commenting step by step because something is wrong
    f(t) = RGB(t, t, t) # lazy grayscale RGB
    image1 = graybrotH(system, window, its) # get preimage 
    transp = !isMoreHorizontal(image1)
    if transp
        image1 = transpose(image1)
    end  
    image = f.(image1) #color it in gray scale 
    inds = getPointsOfInterest(points, image1) # get however many points 
    z(i, j) = indexToComplex(window, i, j) # what it says 

    interpolator = [ inds[i][3-j] for i = 1:points, j = 1:2] # set up array to interpolate
    # interpolating by index 
    vec = polyRegress(interpolator, degree) # poly coeffs from regression 
    
    f = polyOffset( vec , 0) # polynomial IN INDEX

    (yRes, xRes) = size(image) # get res 
    xCent = real(window.center) # complex value of center 
    dx = window.width/xRes # pixel width 
    xInd(x) = (x-xCent)/dx + xRes/2 
                    # divide difference from center by dx for pixels, then add to center pixel

    # take x, get x index. apply to interpolation,
    # plot as point in complex plane 
    func(x) = imag(z(f(xInd(x)),xInd(x))) 
    

    yCent = imag(window.center)
    dx = window.width/xRes
    yInd(y) = floor(Int, div((yCent-y),dx) + div(yRes, 2)) #yindex function in same way

    left = xCent - window.width/2 #left edge of screen 
    x(i) = left + (i-1)*dx #takes index, returns x value 
    lastone = yInd(func(x(1)))
    for i = 1:xRes
        yI = yInd(func(x(i)))
        if 1 <= yI && yI <= yRes
            if lastone < yI
                for l = lastone:yI
                    if (0 < l) && (l<=yRes)
                        image[l, i] = RGB(0,0.8,1)
                    end
                end
            else 
                for l = yI:lastone
                    if (0 < l) && (l<=yRes)
                        image[l, i] = RGB(0,0.8,1)
                    end
                end
            end 
        end 
        lastone = yI
    end 
    for (i, j) in inds
        for k = (i-2):(i+2)
            for l = (j-2):(j+2)
                if (0 < k) && (k<=yRes) &&(0 < l) && (l<=xRes)
                    image[k, l] = RGB(1, 0.5, 0)
                end
            end
        end 
    end
    transp ? transpose(image) : image 
end 



function julianimate(frames, points, degree, system, window, its)
    f = interpolatePointsOfInterest(points, degree, system, window, its)
    jWindow = FractalWindow(300, 300, 0, 0.4)
    width = window.width
    center = window.center
    xStart = real(center) - width/2
    x(i) = xStart + i*width/frames 
    z(i) = x(i) + f(x(i))*im
    for i in 1:frames 
        jul = render(julia(z(i)), jWindow, 600, sBubblegum2)
        save(
            "julianimation\\"*threeDigiter(i)*".png", 
            map(clamp01nan, jul))
        
    end 
end 

function julianimdemo(frames, points, degree, system, window, its)
    f = interpolatePointsOfInterest(points, degree, system, window, its)
    
    mWindow = FractalWindow(600, 600, window.center, window.width)
  
    width = window.width
    center = window.center
    xStart = real(center) - width/2
    # commenting step by step because something is wrong
    fc(t) = RGB(t, t, t) # lazy grayscale RGB
    image1 = graybrotH(system, mWindow, its) # get preimage 

    image = fc.(image1) #color it in gray scale 
    inds = getPointsOfInterest(points, image1) # get however many points 
    z(i, j) = indexToComplex(mWindow, i, j) # what it says 
  
    interpolator = [ inds[i][3-j] for i = 1:points, j = 1:2] # set up array to interpolate
      # interpolating by index 
    vec = polyRegress(interpolator, degree) # poly coeffs from regression 
      
    f = polyOffset( vec , 0) # polynomial IN INDEX
  
    (yRes, xRes) = size(image) # get res 
    xCent = real(mWindow.center) # complex value of center 
    dx = window.width/xRes # pixel width 
    xInd(x) = (x-xCent)/dx + xRes/2 
                      # divide difference from center by dx for pixels, then add to center pixel
  
    # take x, get x index. apply to interpolation,
    # plot as point in complex plane 
    func(x) = imag(z(f(xInd(x)),xInd(x))) 
      
      
    yCent = imag(mWindow.center)
    dx = mWindow.width/xRes
    yInd(y) = floor(Int, div((yCent-y),dx) + div(yRes, 2)) #yindex function in same way
  

    jWindow = FractalWindow(600, 600, 0, 0.4)
    width = window.width
    center = window.center
    xStart = real(center) - width/2
    x(i) = xStart + i*width/frames 
    z(i) = x(i) + func(x(i))*im
    
    lastone = yInd(func(x(1)))
    for (i1, j) in inds
        for k = (i1-2):(i1+2)
            for l = (j-2):(j+2)
                if (0 < k) && (k<=yRes) &&(0 < l) && (l<=xRes)
                    image[k, l] = RGB(1, 0.5, 0)
                end
            end
        end 
    end
    for i in 1:frames 
        jul = render(julia(z(i)), jWindow, 600, sBubblegum2)

        yI = yInd(func(x(i)))
        if 1 <= yI && yI <= yRes
            if lastone < yI
                for l = lastone:yI
                    if (0 < l) && (l<=yRes)
                        image[l, i] = RGB(0,0.8,1)
                    end
                end
            else 
                for l = yI:lastone
                    if (0 < l) && (l<=yRes)
                        image[l, i] = RGB(0,0.8,1)
                    end
                end
            end 
        end  
        lastone = yI
        out = hcat(image, jul)
        save(
            "julianimationdemo\\"*threeDigiter(i)*".png", 
            map(clamp01nan, out))
    end 
end 
begin 
fullMWindow = FractalWindow(200, 200, -1.765 +0.01im, 0.1)
z = randZoomPoint(Mandelbrot, 2, fullMWindow, 5, 10000)
newWindow = FractalWindow(300, 300, z, 0.005)
a = render(Mandelbrot, newWindow, 10000)
end

julianimdemo(600, 50, 4, Mandelbrot, newWindow, 1000)