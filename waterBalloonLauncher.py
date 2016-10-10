GlowScript 2.1 VPython
"""
Launches a water balloon until it hits the ground, using 2nd-order explicit method
Compares the effect of different approximations/assumptions:
*No drag, mass of launcher negligible, no gravity in the launcher

Persistent assumptions: 
*Elastic force is always at angle = ang ( = 45deg)
*Water balloon is a sphere
*Ignores density and thickness of the water balloon rubber.
*Stops when it gets to its original height, at the beginning in the launcher!!

The rate doesn't affect printed impact point
"""

#To do: 
#graphing posx/y, velx/y, accx/y, Eg/k/tot
#
Li = 1.25     #arange(L0+0.1,2.0,0.02), look ~line 33 for
m = 0.175         #0.15 to 0.2 kg
ang = radians(45)        #degrees, gets converted to radians

dt = 2**-7
the_rate=500

makeGraphs = True
#These are cycled later
#ignoreDrag = True   #ignore drag: True-->goes farther
#ignoreGinL = True   #ignore force of gravity in the launcher: True-->goes farther
#ignoreML = True    #ignore work done by bands on launcher (mL): True-->goes farther

canv_width=1500

#Launcher properties
#elastic constant of launcher, measured
k = 91.68       #Measured versus band length, F = k(L - L0) DOUBLE CHECK
mL = 0.188      #kg
cmL = 0.523     #center of mass, factor of launcher
L0 = 0.7173     #m, unstretched length of single band
m_ratio = m/(m+mL*(cmL**2))   #based off of work split between m and mL

#Initial water balloon height
wby0 = 0        #meters above ground, where it starts in the launcher

#Ground properties
xmax = 25       #m
gthick = 1

scene=display(width=canv_width,height=(9/16)*canv_width,
    center=vec(0.5*xmax,wby0,0),autoscale=True,align="left")

ground = box(pos=vec(0.5*xmax,-0.5*gthick,0), length=xmax,width=10,
            height=gthick,color=color.green)

#Water Ball (wb) properties, including drag
rho_w = 1000                #kg/m^3, mass density of water
vol = m/rho_w               #m^3, volume of water balloon
R = (3*vol/(4*pi))**(1/3)   #radius of water balloon
A = pi*R**2                 #cross sectional area
Cd = 0.5                    #depends on shape and Reynold's number, NASA has val
rho_a = 1                   #kg/m^3, mass density of air
drag_factor = 0.5*rho_a*A*Cd

sq2 = 0.5*sqrt(2)
unitVecAng = vec(cos(ang),sin(ang),0)
the_colors = [color.orange,color.green,color.yellow,color.blue,color.red]
ignore_mask = [[True,True,True],
                [False,True,True],
                [True,False,True],
                [True,True,False],
                [False,False,False]]
[ignoreDrag,ignoreGinL,ignoreML] = ignore_mask[0]
wb=sphere(pos=vec(0,wby0,0), radius=R, color=the_colors[0], make_trail=True)
wb.m = m
#Other constants
g = 9.81        #m/s^2
Fg = vec(0,-m*g,0)


#Create graphs
if makeGraphs:
    n_plot=1
    posx_g = graph(width=550,height=300,align="left")
    posy_g = graph(width=550,height=300,align="right")
    velx_g = graph(width=550,height=300,align="left")
    vely_g = graph(width=550,height=300,align="right")
    accx_g = graph(width=550,height=300,align="left")
    accy_g = graph(width=550,height=300,align="right")

#z = Straight line stretch (path water balloon follows), meters
#s = band length, the useful thing to measure, meters
z = sqrt(Li**2 - L0**2)
yrelease = z*sin(ang)
xrelease = z*cos(ang)
#attachment points for thw launcher band (want to simulate this for effect)
p1 = vec(xrelease,yrelease,-L0)
p2 = vec(xrelease,yrelease,L0)
band1=cylinder(pos=p1,axis=wb.pos-p1,radius=R,color=color.white)
band2=cylinder(pos=p2,axis=wb.pos-p2,radius=R,color=color.white)

for i in range(5):
    t = 0
    n = 0
    wb.v = vec(0,0,0)           #released from rest   
    #reset waterballoon for next launch, different approximations
    wb.make_trail = False
    wb.pos = vec(0,wby0,0)
    #Change color and approximations
    wb.color = the_colors[i]
    wb.trail_color=wb.color
    wb.make_trail = True
    [ignoreDrag,ignoreGinL,ignoreML]=ignore_mask[i]
    if makeGraphs:
        posx_s = series(graph=posx_g,label='x (m)',color=the_colors[i])
        posy_s = series(graph=posy_g,label='y (m)',color=the_colors[i])
        velx_s = series(graph=velx_g,label='v_x (m/s)',color=the_colors[i])
        vely_s = series(graph=vely_g,label='v_y (m/s)',color=the_colors[i])
        accx_s = series(graph=accx_g,label='a_x (m/s^2)',color=the_colors[i])
        accy_s = series(graph=accy_g,label='a_y (m/s^2)',color=the_colors[i])
    while wb.pos.y > wby0 - R:      #stop @ original height, at beginning of launch, NOTE!
        rate(the_rate)
        
        Ftotal = vec(0,0,0)
        Fdrag = (-drag_factor*mag(wb.v))*wb.v
        
        #wb starts at y = 0, moves up until released
        if wb.pos.y < yrelease and wb.pos.x < xrelease:
            #z2 = mag2(wb.pos - vec(yrelease,yrelease,0))
            #L = sqrt(z2 + L0**2)
            band1.axis=wb.pos-p1
            band2.axis=wb.pos-p2
            L=mag(band1.axis)
            FL = (k*(L - L0))*unitVecAng
            if not ignoreML:
                FL = FL*m_ratio  #water balloon's mass-fraction of system
            
            Ftotal = Ftotal + FL
            if not ignoreGinL:
                Ftotal = Ftotal + Fg
            if not ignoreDrag:
                Ftotal = Ftotal + Fdrag
        else:
            Ftotal = Ftotal + Fg
            if not ignoreDrag:
                Ftotal = Ftotal + Fdrag
        wb.a = Ftotal/wb.m
        
        #Graph
        if makeGraphs and n%n_plot == 0:
            #Plot that it gets released at x = 0
            posx_s.plot(t,wb.pos.x-xrelease)
            posy_s.plot(t,wb.pos.y)
            velx_s.plot(t,wb.v.x)
            vely_s.plot(t,wb.v.y)
            accx_s.plot(t,wb.a.x)
            accy_s.plot(t,wb.a.y)
        
        #update using midpoint (central point approximation to derivative)
        if n == 0:  #first step, use euler
            vold = wb.v
            rold = wb.pos
            dt *= 0.5
            
        rtmp = rold + wb.v*2*dt
        vtmp = vold + wb.a*2*dt
        
        #store old point for next iteration
        rold = wb.pos
        vold = wb.v
        #update position and velocity
        wb.pos = rtmp
        wb.v = vtmp
        
        if n == 0:
            dt *= 2
        
        n += 1
        t += dt
    
    print('Traveled',wb.pos.x-xrelease,
            '(m) from',Li,'(m) band-length while ignoring:',
            '[drag,gravity in launcher, mass of launcher] =',
            [ignoreDrag,ignoreGinL,ignoreML])
