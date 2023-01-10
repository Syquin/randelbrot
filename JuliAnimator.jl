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
        imfilter(img, Kernel.Gaussian(10)),
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
    image = prerender(system, window, its)
    inds = getPointsOfInterest(points, image)
    z((i, j)) = indexToComplex(window, i, j)
    parts(comp, i) = (i == 1) ? real(comp) : imag(comp)
    interpolator = [ parts(z(inds[i]),j) for i = 1:points, j = 1:2]
    vec = polyRegress(interpolator, degree)
    println(vec)
    polyOffset( vec , 0)
end 

function demonstrateInterpolate(points, degree, system, window, its)
    f(t) = RGB(t, t, t)
    image1 = prerender(system, window, its)
    image = f.(image1)
    inds = getPointsOfInterest(points, image1)
    z((i, j)) = indexToComplex(window, i, j)
    parts(comp, i) = (i == 1) ? real(comp) : imag(comp)
    println(z.(inds))
    interpolator = [ parts(z(inds[i]),j) for i = 1:points, j = 1:2]
    vec = polyRegress(interpolator, degree)
    func = polyOffset( vec , 0)
    (yRes, xRes) = size(image)
    xCent = real(window.center)
    yCent = imag(window.center)
    dx = window.width/xRes
    yInd(y) = floor(Int, div((yCent-y),dx) + div(yRes, 2))
    left = xCent - window.width/2 
    x(i) = left + dx*(i-1)
    for i = 1:xRes
        x1 = x(i)
        y1 = func(x1)
        yI = yInd(y1)
        if 1 <= yI && yI <= yRes
            image[yI, i] = RGB(1,0,0)
        end 
    end 
    for (i, j) in inds
        for k = (i-2):(i+2)
            for l = (j-2):(j+2)
                if (0 < k) && (k<=yRes) &&(0 < l) && (l<=xRes)
                    image[k, l] = RGB(0, 1, 0)
                end
            end
        end 
    end
    image 
end 



function julianimate(frames, points, degree, system, window, its)
    f = interpolatePointsOfInterest(points, degree, system, window, its)
    jWindow = FractalWindow(300, 300, 0, 3)
    width = window.width
    center = window.center
    xStart = real(center) - width/2
    x(i) = xStart + i*width/frames 
    z(i) = x(i) + f(x(i))*im
    for i in 1:frames 
        jul = render(julia(z(i)), jWindow, 600)
        @time save(
            "julianimation\\"*threeDigiter(i)*".png", 
            map(clamp01nan, jul))
    end 
end 

fullMWindow = FractalWindow(200, 200, -0.3 + 0.65im, 0.1)
render(Mandelbrot, fullMWindow, 10000)

julianimate(300, 30, 5, Mandelbrot, fullMWindow, 1000)