using Images
using ImageIO

# take a list of numbers and turns them into the corresponding polynomial
function polyOffset(coefs::Vector, offset::Real)
    c = reverse(coefs)
    terms = x ->  [(c[i]*((x-offset)^(i-1))) for i = 1:length(c)]
    return sum ∘ terms
end 

# take a list of points xy expressed as a n by 2 matrix and finds the polynomial of degree n of best fit
function polyRegress(xy, n)
    A = ones(typeof(0.0), (length(xy[:,1]), n+1))
    i=1
    for v in xy[:,1]
        x = v[1]
        A[i, 1:n] = [x^(n-t+1) for t = 1:n]
        i += 1
    end
    At = transpose(A)
    b = xy[:, 2]
    return inv(At*A)*At*b
end


# takes RGB values and creates a color gradient from black to your color (r,g,b) 
# to white using quadratic interpolation. 
function b2c2w(r, g, b)  
    test1 = [0 0; 0.3 r; 1 1]
    test2 = [0 0; 0.3 g; 1 1]
    test3 = [0 0; 0.3 b; 1 1]

    testParaV1 = polyRegress(test1, 2)
    testParaV2 = polyRegress(test2, 2)
    testParaV3 = polyRegress(test3, 2)

    testPara1 = polyOffset(testParaV1, 0)
    testPara2 = polyOffset(testParaV2, 0)
    testPara3 = polyOffset(testParaV3, 0)


    testGrad(t) = RGB(testPara1(t),testPara2(t),testPara3(t) )
end 


#  COLOR GRADIENTS BEGIN
# Color gradients are expressed as functions from the unit interval [0,1] to color spPurple

wblue = b2c2w(0.43, 0.25, 0.36)
wpink = b2c2w(0.16, 0.4, 0.46)
sWildberry(t) =  t == 0 ? RGB(0.0,0.0,0.0) : ( t < 0.5 ? wblue((2t)) : wpink((2-2t)) )


tpink = b2c2w(0.63, 0.4, 0.5)
tblue = b2c2w(0.16, 0.45, 0.66)
sBubblegum2(t) =  t == 0 ? RGB(0.0,0.0,0.0) : ( t < 0.5 ? tblue((2t)) : tpink((2-2t)) )
sBubblegum2.(0:0.01:1)


nbpurp = b2c2w(0.4, 0, 0.6)
nbyellow = b2c2w(0.5, 0.5, 0)
sNB(t) =  t == 0 ? RGB(0.0,0.0,0.0) : ( t < 0.5 ? nbpurp((2t)) : nbyellow((2-2t)) )
sNB.(0:0.01:1)


apurp = b2c2w(0.4, 0, 0.6)
agray = b2c2w(0.1, 0.1, 0.1)
sPurple(t) =  t == 0 ? RGB(0.0,0.0,0.0) : ( t < 0.5 ? apurp((2t)) : agray((2-2t)) )
sPurple.(0:0.01:1)

agreen = b2c2w(0.25, 0.5, 0.25)
sGreen(t) =  t == 0 ? RGB(0.0,0.0,0.0) : ( t < 0.5 ? agreen((2t)) : agray((2-2t)) )

orange = b2c2w(0.7, 0.25, 0)
orange.(0:0.01:1)
pink = b2c2w(0.43, 0.08, 0.5)
pink.(0:0.01:1)
sRed(t) =  t == 0 ? RGB(0.0,0.0,0.0) : ( t < 0.5 ? orange((2t)) : pink((2-2t)) )
sRed.(0:0.01:1)


green = b2c2w(0.2, 0.5, 0.3)
blue = b2c2w(0.23, 0.1, 0.66)
sBlue(t) =  t == 0 ? RGB(0.0,0.0,0.0) : ( t < 0.5 ? blue((2t)) : green((2-2t)) )
sBlue.(0:0.01:1)


sThermal(t::Float64) = t == 0 ? RGB(0.0,0.0,0.0) : RGB(t, t^2, 3*t*(1-t)^2)

sNuclear(t::Float64) = t == 0 ? RGB(0.0,0.0,0.0) : RGB(t^2, t*0.9, (2.6*sqrt(t)*(1-sqrt(t)))^(4))

sBrass(t::Float64) = t == 0 ? RGB(0.0,0.0,0.0) : RGB(t^1.4, t*0.9, (2.6*sqrt(t)*(1-sqrt(t)))^(3))


deepSky = RGB(0.1, 0.125, 0.5)
rust = RGB(0.7, 0.4, 0)

convexCombo(x, y, t) = (1-t)*x + t*y

convexComboS(x, y, t) = x*(1-t) + y*t

fadeToWhite(x, y, z) = t -> RGB(convexComboS(x, 1, t), convexComboS(y, 1, t), convexComboS(z, 1, t) )

fadeToWhite(0.9, 0.6, 0).(0:0.01:1)


blueGrad = b2c2w(0.1, 0.125, 0.5)
orangeGrad = b2c2w(0.8, 0.6, 0)

sWiki(t) =  t == 0 ? RGB(0.0,0.0,0.0) : ( t < 0.5 ? blueGrad((2t)) : orangeGrad((2-2t)) )
sWiki.(0:0.01:1)


sBY(t::Float64) = t == 0 ? RGB(0.0,0.0,0.0) : RGB(t^1.4, t*0.9, (3*sqrt(t)*(1-sqrt(t)))^(3))
sBY.(0:0.001:1)
spikeHX(t::Float64) = if (t <= 0.45) 0 elseif t <= 1/2 400*(t-0.45)^2 elseif t <= 0.55 400*(0.05- (t-0.5))^2 else 0  end

sEyeCandy(t::Float64) = t == 0 ? RGB(0.0,0.0,0.0) : RGB(4*t*(1-t),  1-t , (1-t))

sEyeCandy.(0:0.001:1)

sBubblegum(t::Float64) = t == 0 ? RGB(0.0,0.0,0.0) : RGB(sqrt(3*t*(1-t)), t, t^(0.75))
sBubblegum.(0:0.01:1)

spikeH(t::Float64) = if (t <= 1/3) 0 elseif (t <= 1/2) sqrt(6*(t-1/3)) elseif (t <= 5/6) 1 else sqrt(6-6t) end
spike(t::Float64) = float(spikeH(t % 1))


sRainbow(t::Float64) = t==0 ? RGB(0.0,0.0,0.0) : RGB(spike(t+1/3),spike(t),spike(t-1/3))
sRainbow.(0:0.001:1)
