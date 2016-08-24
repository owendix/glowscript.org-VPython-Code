# -*- coding: utf-8 -*-
GlowScript 2.1 VPython
#Lesson 9 Impulse momentum

scene.align="left"

nVecDraw = 700
makeGraphs = True
#Create graphs
if makeGraphs:
    n_plot=200
    posx_g = graph(width=500,height=300,align="left")
    posy_g = graph(width=500,height=300,align="right")
    velx_g = graph(width=500,height=300,align="left")
    vely_g = graph(width=500,height=300,align="right")
    accx_g = graph(width=500,height=300,align="left")
    accy_g = graph(width=500,height=300,align="right")
    posx_s = series(graph=posx_g,label='x (m)',color=color.red)
    posy_s = series(graph=posy_g,label='y (m)',color=color.red)
    velx_s = series(graph=velx_g,label='v_x (m/s)',color=color.blue)
    vely_s = series(graph=vely_g,label='v_y (m/s)',color=color.blue)
    accx_s = series(graph=accx_g,label='a_x (m/s^2)',color=color.black)
    accy_s = series(graph=accy_g,label='a_y (m/s^2)',color=color.black)


def drawmotionvecs(obj,origin=vec(0,-1.5,0)):
    #These objects have methods: pos, p, F, m
    #So far this seems easier than
    #the build in make_arrow feature, which has little
    #documentation
    C = 0.002    #scale for shaftwidth
    rmag = mag(obj.pos)
    sw = C*rmag
    #rpscale = 1
    rvec = arrow(pos=origin,axis=obj.pos-origin,shaftwidth=sw,color=color.red)
    pvec = arrow(pos=obj.pos,axis=obj.p/obj.m,shaftwidth=5*sw,color=color.blue)
    Fvec = arrow(pos=obj.pos,axis=obj.F/obj.m,shaftwidth=10*sw,color=color.white)
    #Momentum has too widely varying scale compared to r to map both and 
    #have significance in the momentum magnitude
    return [rvec,pvec,Fvec]
    
    #return rvec

scene = display(forward=vec(0,-0.5,-1))
scene.fullscreen = True

ball = sphere(pos=vec(-60,2,0),radius=4,color=vec(1,0,0),
    make_trail=True)#,trail_type="points",interval=200)
ground = box(pos=vec(0,-2,0), length=120,height=1,width=20,color=vec(0,1,0))

ball.m = 5 #kg
Fnet = vec(0,-50,0) #N
ball.p = vec(120,120,0)
ball.F = Fnet    #can multiply a scalar (m^-1) by a vector in VPython

#balltrail = curve(color=ball.color)
#Try to replace with built in feature

t=0
dt=0.001

    
motionvecs=[]
n = 0

while ball.pos.x < 60:
    rate(5000)

    if n%nVecDraw==0:
        bvecs = drawmotionvecs(ball)
        #for b in bvecs:
        #    b.color=ball.color
        motionvecs.extend(bvecs)
    #Graph
    if makeGraphs:
        if n%n_plot == 0:
            posx_s.plot(t,ball.pos.x)
            posy_s.plot(t,ball.pos.y)
            v=ball.p/ball.m
            velx_s.plot(t,v.x)
            vely_s.plot(t,v.y)
            a=ball.F/ball.m
            accx_s.plot(t,a.x)
            accy_s.plot(t,a.y)
        
    ball.p = ball.p + ball.F*dt
    
    ball.pos=ball.pos+ball.p*(dt/ball.m)
    
    #balltrail.append(ball.pos)
    n += 1
    t += dt