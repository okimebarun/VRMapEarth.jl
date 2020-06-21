#=
title: VRearth( manipulatable virtual earth map ) 
version: 0.1.0
date: 2020/06/21
author: okimebarun
url: https://github.com/okimebarun/
url: https://qiita.com/oki_mebarun/
=#

module VRMapEarth

greet() = print("Hello World!")

############################################################################
using Interact
using Plots
using Images, ImageView, ImageMagick
using Blink

#=
 calculation part
=#

# parameters 
mutable struct ParamDict
    p::Float32  # z-axis of eye point
    r::Float32  # sphere radius
    v0::Int     # mid-point of v
    w0::Int     # mid-point of w
    m0::Int     # mid-point of m
    n0::Int     # mid-point of n
    z1::Float32 # z-axis of img-plane
    z2::Float32 # z-axis of background-plane (dont use now)
    msh::Int    # map size for longtitude
    msv::Int    # map size for latitude
end;

# mapping functions

# get (a,b) from (v,w)
function imcoef(d::ParamDict, v::Int, w::Int)
    t=d.p - d.z1
    a=Float32(v - d.v0)/t
    b=Float32(w - d.w0)/t
    (a,b)
end

# get (flag, x,u,z) from (a,b)
function spxyz(d::ParamDict, a::Float32, b::Float32)
    ab1 = a^2+b^2+1
    r1 = d.p^2 -(d.p^2 - d.r^2)*ab1
    if r1 < 1e-10
        (false, 0.0, 0.0, 0.0)
    else
        t = (d.p - sqrt(r1)) / ab1
        x = a*t
        y = b*t
        z = -t + d.p
        (true, x, y, z)
    end
end

# get (long, lat ,coef) from (x,y,z) and (long0, lat0)
function polerax(d::ParamDict, x::Float32, y::Float32, z::Float32,
        long0::Int, lat0::Int)
    
    phi = pi*(Float32(lat0)/180.0+1.0)
    ax = 0
    ay = cos(phi)
    az = sin(phi)
    polerv = [ax, ay,  az] # earth's axis
    polerh = [ 0, az, -ay] # orthogonal to polerv
    
    # get latitude
    r = d.r
    c1 = [x,y,z]'*polerv/r
    lat = acos(c1*0.9999)

    # get longtitude
    # inline function
    cross3d(a,b)=[a[2]b[3]-a[3]b[2], a[3]b[1]-a[1]b[3], a[1]b[2]-a[2]b[1]]
    b  = [x,y,z]/r - polerv*c1 # projection on polerv
    b  = b/sqrt(b'*b) # unitilize
    cv = cross3d(b, polerh) # outer product
    c2 = cv'*polerv
    c3 = b' *polerh
    long = acos(c3*0.9999) * sign(c2)+ pi*Float32(long0)/180.0   
    
    # convert to 360 / 180 degree
    long2 = Float32(360.0*(long/(2pi) + 2.0) % 360.0) # Mod on float
    lat2  = Float32(180.0*( lat/(pi)  + 2.0) % 180.0)

    # reflection coefficient
    coef = [x, y, z]'*[0.0, 0.0, 1.0]/d.r
    (long2, lat2, coef)
end

# get (mx,my) from (long, lat)
function mapxy(d::ParamDict, long::Float32, lat::Float32)
    mx = round(Int, long/360.0*Float32(d.msh-1))+1
    my = round(Int, lat /180.0*Float32(d.msv-1))+1
    (mx, my)
end

# get (m, n) from (v,w)
function backax(d::ParamDict, v::Int, w::Int)
    m = round(Int, ((v-d.v0)*d.m0)/d.v0 + d.m0)+1
    n = round(Int, ((w-d.w0)*d.n0)/d.w0 + d.n0)+1
    (m, n)
end

# get (flag, my, mx, coef) from (w, v, long0, lat0)
function convToMap(d::ParamDict, w::Int, v::Int, long0::Int, lat0::Int)
    (a,b) = imcoef(d, v, w)
    (bon, x, y, z) = spxyz(d, a, b)
    if bon
        (long, lat, coef) = polerax(d, x,y,z, long0, lat0)
        (mx, my) = mapxy(d, long, lat)
        (true, my, mx, coef)
    else
        (m, n) = backax(d, v, w)
        (false, n,  m, 0.0)
    end
end

#=
 GUI part
=#

# main
function real_body()
    # load images 
    img_galaxy = load("img_galaxy.jpg")
    img_earth  = load("img_earth.jpg")
    
    # setup parameter
    d1 = ParamDict(5000.0, 300.0, 400, 300, 500, 300, 510.0, -600.0, 1000, 800)
    d1.msh = size(img_earth)[2]
    d1.msv = size(img_earth)[1]
    d1.m0 = round(Int,size(img_galaxy)[2]/2-1)
    d1.n0 = round(Int,size(img_galaxy)[1]/2-1)
    
    # view plane
    vs=1:800
    ws=1:600
    
    # functions
    # convert to image
    function imconv1(p)
        # shade effect
        rgb0 = RGB{N0f8}(0.1,0.1,0.2)
        # select pixel from two images depending on the first flat
        p[1] ? RGB{N0f8}( img_earth[p[2],p[3]]*p[4] + rgb0*(1.0-p[4]) ) : img_galaxy[p[2],p[3]]
    end
    
    # visualize the result
    function view1(long::Int, lat::Int, z1::Float32, d1=d1, ws=ws, vs=vs)
        d1.z1 = z1
        imap1a(w,v)=convToMap(d1, w, v, long, lat)
        # Julia like mapping!
        m   = imap1a.( ws, vs' )
        img = imconv1.(m)
        colorview(RGB, img)
    end
    
    # GUI objects
    long = slider(0:10:360,  label="longitude");
    lat  = slider(-90:5:90,  label="latitude ");
    zoom = slider(-5000:100:4500,label="zoom in/out");
    slider01 = vbox(long,lat,zoom);
    
    # interacting function
    img_plot01=map(
        (long, lat, zoom) ->view1(long, lat, Float32(zoom)),
        map(observe, [long, lat, zoom])...,
    );
    
    # UI
    ui = vbox(slider01, img_plot01);
    
    # show window via Blink
    w = Window(Dict(:title=>"VR earth"));
    size(w,818,800);
    body!(w, ui);
    
    while( active(w))
        sleep(1)
    end
end
#

############################################################################

#
function julia_main()
    try
        real_main()
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end

#
function real_main()
    @show ARGS
    @show Base.PROGRAM_FILE

    # do something
    real_body()
	#
	
    return
end

if abspath(PROGRAM_FILE) == @__FILE__
    real_main()
end
#

end # module
