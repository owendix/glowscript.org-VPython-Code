# -*- coding: utf-8 -*-
GlowScript 2.1 VPython
scene = display(forward=vec(0,-0.5,-1),fullscreen=True,align="left")

nVecDraw = 1000
makeGraphs = True
#Create graphs
if makeGraphs:
    n_plot=200
    posx_g = graph(width=550,height=300,align="left")
    posy_g = graph(width=550,height=300,align="right")
    velx_g = graph(width=550,height=300,align="left")
    vely_g = graph(width=550,height=300,align="right")
    accx_g = graph(width=550,height=300,align="left")
    accy_g = graph(width=550,height=300,align="right")
    posx_s = series(graph=posx_g,label='x (m)',color=color.red)
    posy_s = series(graph=posy_g,label='y (m)',color=color.red)
    velx_s = series(graph=velx_g,label='v_x (m/s)',color=color.blue)
    vely_s = series(graph=vely_g,label='v_y (m/s)',color=color.blue)
    accx_s = series(graph=accx_g,label='a_x (m/s^2)',color=color.black)
    accy_s = series(graph=accy_g,label='a_y (m/s^2)',color=color.black)

def drawmotionvecs(obj,origin=vec(0,-19.5,0)):
    #These objects have methods: pos, v, a
    #So far this seems easier than
    #the build in make_arrow feature, which has little
    #documentation
    C = 0.002    #scale for shaftwidth
    rmag = mag(obj.pos)
    sw = C*rmag
    #rpscale = 1
    rvec = arrow(pos=origin,axis=obj.pos-origin,shaftwidth=sw,color=color.red)
    vvec = arrow(pos=obj.pos,axis=obj.v,shaftwidth=5*sw,color=color.blue)
    avec = arrow(pos=obj.pos,axis=obj.a,shaftwidth=10*sw,color=color.white)
    #Momentum has too widely varying scale compared to r to map both and 
    #have significance in the momentum magnitude
    return [rvec,vvec,avec]
    
    #return rvec

ball = sphere(pos=vec(-70,-16,0),radius=4,color=vec(1,0,0),
    make_trail=True)#,trail_type="points",interval=300)
ground = box(pos=vec(0,-20,0), length=140,height=1,width=20,color=vec(0,1,0))

ball.v = vec(18,35,0)
ball.a = vec(0,-9.8,0)

#balltrail = curve(color=ball.color)
#Try to replace with built in feature

t=0
dt=0.001

motionvecs=[]
n = 0
while ball.pos.y > -18:
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
            velx_s.plot(t,ball.v.x)
            vely_s.plot(t,ball.v.y)
            accx_s.plot(t,ball.a.x)
            accy_s.plot(t,ball.a.y)
    
    ball.v = ball.v + ball.a*dt
    
    ball.pos = ball.pos + ball.v*dt
    
    #balltrail.append(ball.pos)
    n += 1
    t += dt