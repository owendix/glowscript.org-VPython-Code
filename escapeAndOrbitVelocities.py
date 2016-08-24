# -*- coding: utf-8 -*-
GlowScript 2.1 VPython

#Earth properties
solidEarth = False   #if True object cannot fall through the Earth
                    #False demonstrates that the object would be an ellipse
                    #but only if earthStructure = ignore or const_g
uniformEarth = True #True assumes a uniform density Earth, False uses g=GM/r^2
"""
Fun side note:
In 2015, Alexander Klotz showed the uniform Earth gives less accurate diameter
travel time through the center of the Earth than when you consider a constant 
9.81m/s^2 toward the center, given the layered structure of the Earth and its 
corresponding density:
http://scitation.aip.org/content/aapt/journal/ajp/83/3/10.1119/1.4898780
The more accurate answer is ~38min, which matches a constant g=9.81m/s^2
The uniform Earth yields a differential equation:
dr/dt = -(GM/R^3)r, which has solution (starting from r = R, v = 0) or
r(t) = Rcos(w*t), where w = sqrt(GM/R^3)
Since w = 2*pi/T, and T is the period, this yields a period of 
T = 2*pi*sqrt(R^3/GM) = 2*pi*sqrt(R/9.81m/s^2) = 2*pi*(6.371e6m/9.81m/s^2)
T = 84min --> going a full diameter is T/2 = 42min
*Uniform gravity, g = 9.81m/s^2 
--> You're actually traveling half at -g, half at positive g and this
trip is symmetric so we'll find the time it takes to go R then double that
time: R = (gt^2)/2
t_full2R = 2*sqrt(2*R/g) = 2*sqrt(2*6.371e6m/9.81m/s^2) = 37.99min
The two answers are actually off by ~11seconds but the const_g=9.8m/s^2 is closest
*The assumption of a constant g makes three wrong assumptions that cancel pretty closely:
It ignores the fact that as you get closer to the center of the Earth, you are 
being harder because your distance is closer to the Earth's center of mass. But 
it also ignores the fact that the mass of the shell outside of your current position,
as you fall, doesn't contribute to the pull of gravity on you, and the fact that 
there is more dense matter pulling you from closer to the core. Uniform gravity 
correctly ignores the outer shell and the closer distance to the center of mass right, 
but the constant-g's model ignoring the greater Earth density closer to the core 
apparently cancels with the other two assumptions, given the Earth's specific structure.
"""

dt = 10
t_stop = dt*2000
the_rate=100
n_plot = 10
n_stuck = 100

scene.align="left"

#Initial angle, in radians (from horizontal right)
ang = 0
#Initial speed, m/s
v0s = [4000,5000,6000,7000,7025,7050,7075,8000,10000,12000]

#how much smaller to display the ball
rbfact = 0.1       


if solidEarth:
    opac=1
else:
    opac=0.35

M = 5.972e24        #kg
G = 6.67408e-11     #Nm^2/kg^2
GM = G*M
Re = 6.371e6        #m


E_gph=graph(width=450,align="right")
Egs = series(graph=E_gph,color=color.green,label='Eg (J)')
Eks = series(graph=E_gph,color=color.blue,label='Ek (J)')
Etots = series(graph=E_gph,color=color.red,label='Etot (J)')

earth = sphere(pos=vec(0,0,0),radius=6.371e6, color=vec(0.4,0.4,1),opacity=opac)

#ball
m = 1               #kg
m2 = 0.5*m
GMm = GM*m          
rb = Re*rbfact      #m
r_solid = Re*(1 + rbfact)    #s
r_start = r_solid + rb

if solidEarth:
    ball = sphere(pos=vec(0,r_start,0),radius=rb,color=color.red)
else:
    ball = sphere(pos=vec(0,r_start,0),radius=rb,color=color.red,
        make_trail=True)
    
for v in v0s:
    i_stuck = 0
    Egs.delete()
    Eks.delete()
    Etots.delete()
    
    print('v0 =',v,'m/s ...')
    if not solidEarth:
        ball.make_trail=False
        ball.pos = vec(0,r_start,0)
        ball.make_trail=True
    else:
        ball.pos = vec(0,r_start,0)

    ball.v = v*vec(cos(ang),sin(ang),0)
    n = 0
    if v > 9000:
        t_stop = dt*100000
    for t in arange(0,t_stop,dt):
        rate(the_rate)
        
        r2 = mag2(ball.pos)
        r = mag(ball.pos)
        
        #g = GM/r2
        if not uniformEarth or r >= Re:
            Eg = -GMm/r
        elif uniformEarth and r < Re:
            Eg = -GMm/(2*Re)*(3 - (r/Re)**2)  #got it
            """
            Eg = m*Vg
            Gradient(Vg) = - g (vec)
            Vg(r) = -Integral(g*dr) from infinity to r (through R)
            """
            
        Ek = m2*mag2(ball.v)
        Etot = Eg+Ek
        if solidEarth:
            if r <= r_solid:
                print('v0 =',v,'m/s: IMPACT at latitude',
                    atan(ball.pos.y/ball.pos.x)*180/pi,'deg')               
                break
        else:
            if r < 0.05*r_solid:
                if i_stuck < n_stuck:
                    i_stuck += 1
                else:
                    print('v0 =',v,'m/s seems stuck inside Earth')
                    break
        
        if n%n_plot == 0:
            Egs.plot(t,Eg)
            Eks.plot(t,Ek)
            Etots.plot(t,Etot)
        if not uniformEarth or r >= Re:
            acc = -GM/r**3
        elif uniformEarth and r < Re:
            acc = -GM/Re**3
        ball.a = acc*ball.pos        #acceleration pts to center of earth
        
        ball.v = ball.v + ball.a*dt
        ball.pos = ball.pos + ball.v*dt
        n += 1

    if t >= (t_stop-dt) and r > r_solid:
        print('v0 =',v,'m/s entered orbit!')