GlowScript 2.1 VPython
"""
Drop three perfectly elastic balls of different sizes.
Does not work elegantly. rate(frequency,move) --> 1/frequency is linear with 
the number of button pushes, which I discovered after some trial, error, and 
hypothesis testing. My Kluge does fix the problem, as precisely as I can tell. 

Asynchronous calls used but restarting leads 
to an increase in rate(the_rate,move). Since I don't know quite how
rate or wait works I don't know how to fix it. I also tried to include wait as an 
argument to a function that contains it, unless it contains it in rate in which 
case it seems to be included automatically (it must still be included 
as an argument in the calling function statement though. This seems 
to not allow an interruption at all. If I can pass wait with some other 
callback function that may work.
"""

ms = [1, 0.5, 0.1]
rs = [0.1, 0.07937, 0.04642]
space = rs[0]/10
colors = [color.yellow, color.cyan, color.magenta] 
y_pos0 = 1.0        #m
g_thick = rs[0]
g_zlim = 10*g_thick
g_xlim = 40*g_thick
balls = [0, 0, 0]

canv_width=1000
canv_height=canv_width*(9/16)
scene=display(width=canv_width,height=canv_height,center=vec(0,y_pos0,0))

scene.title='Drop three perfectly elastic hard balls\n\n'
ground = box(pos=vec(0,-0.5*g_thick,0),axis=vec(0,1,0),length=g_thick,
            width=2*g_zlim,height=2*g_xlim,color=color.green)

auto_scale = False   #follows ball
allow_unstable = False

g = 9.8          #m/s^2, =0 and =9.8 to see effects (minimal)
dt = 1e-4
the_rate = 1e4
stop = False
eps = 1e-10

def check_ground(balls,xlim,zlim):
    dp = 0
    #ground height at 0
    #edges of ground at x = +- xlim, z = +- zlim
    #right now this ignores corners and edges
    for b in balls:
        #too far out
        yg = b.radius
        if b.pos.y <= yg:
            xgp = xlim + b.radius
            xgn = -xgp 
            if b.pos.x >= xgp or b.pos.x <= xgn:
                continue
            zgp = zlim + b.radius
            zgn = -zgp  
            if b.pos.z >= zgp or b.pos.z <= zgn:
                continue
            if b.pos.x <= xlim and b.pos.x >= -xlim:
                if b.pos.z <= zlim and b.pos.z >= -zlim:
                    #shift back inside by too-far amount
                    b.pos = b.pos + vec(0,yg-b.pos.y,0)
                    b.v.y = -b.v.y      #elastic bounce means v in normal dir flips
                    #should have some sort of shift
                    #b.pos = b.pos + vec(0,b.v.y,0)*del_t
                    dp += 2*b.m*abs(b.v.y)
            #check for corners
            #check for edges

    return [balls,dp]
    
def collide_balls(balls):

    for i,a in enumerate(balls[:-1]):
        for j, b in enumerate(balls[i+1:]):
            ab = a.pos - b.pos
            RaRb2 = (a.radius+b.radius)**2
            if mag2(ab) <= RaRb2:
                M = a.m + b.m
                #back up to where they should have contacted
                #assumes gravitational acceleration negligible during impact
                Mr2 = M*mag2(ab) #Mr2 = M*mag2(ab)
                rv = dot(ab,b.v - a.v)
                vba = mag(b.v - a.v)
                rv0 = rv/vba
                b2 = mag2(ab) - rv0**2
                del_t = (sqrt(RaRb2 - b2) - rv0)/vba
                #need something to account for acceleration during del_t
                #maybe make the balls elastic and create an acceleration?
                a.pos = a.pos - a.v*del_t
                b.pos = b.pos - b.v*del_t
                #compute new velocities, acquired through cm calculation
                a.v = a.v + 2*(b.m/Mr2)*rv*ab       #verified
                b.v = b.v - 2*(a.m/Mr2)*rv*ab
                #how far out from impact they should be
                a.pos = a.pos + a.v*del_t    #scoot them out
                b.pos = b.pos + b.v*del_t    #from impact

    return balls
    

def move():
    global balls, stop
    if stop:
        stop = False #restore to original value before exiting move()
        return
    rate(the_rate,move)     #when move() defined, it ignores or increases the_rate?
    #make them move
    for b in balls:
        b.pos = b.pos + b.v*dt
        b.v = b.v + b.a*dt

    #check for ground collisions
    [balls,dp] = check_ground(balls,g_xlim,g_zlim)
    
    #check for interparticle collisions
    balls = collide_balls(balls)
    y_max = 0
    for b in balls:
        if b.pos.y > y_max:
            y_max = b.pos.y
    if y_max < -g_thick:
        return

def Start_Sim():
    global balls #, the_rate
    
    scene.range=2*g_zlim
    scene.autoscale=auto_scale
    xz_comps = vec(0,0,0)
    try:
        for i, b in enumerate(balls):
            b.visible = False
            if i == 0:
                the_y = y_pos0
            else:
                the_y = old_b_pos.y + rs[i-1] + rs[i] + space
            if allow_unstable:
                xz_comps = (rs[2]*eps)*vec.random() #tiny random perturbation
            b.pos=vec(xz_comps.x,the_y,xz_comps.z)
            b.radius=rs[i]
            b.color=colors[i]
            b.visible=True
            b.m = ms[i]
            b.v = vec(0,0,0)
            b.a = vec(0,-g,0)
            old_b_pos = b.pos
    except:
        for i in range(len(rs)):
            if i == 0:
                the_y = y_pos0
            else:
                the_y = old_b_pos.y + rs[i-1] + rs[i] + space
            if allow_unstable:
                xz_comps = (rs[2]*eps)*vec.random() #tiny random perturbation
            balls[i]=sphere(pos=vec(xz_comps.x,the_y,xz_comps.z),radius=rs[i], 
                color=colors[i],visible=True)
            balls[i].m = ms[i]
            balls[i].v = vec(0,0,0)
            balls[i].a = vec(0,-g,0)
            old_b_pos = balls[i].pos
    
    move()

def Restart_Sim(the_radio):
    global stop
    stop = True
    Start_Sim()
    
rs_radio = button(text='Restart',pos=scene.title_anchor,bind=Restart_Sim)

def Auto_Scale(the_radio):    
    global auto_scale, stop
    stop = True
    auto_scale = not auto_scale
    if auto_scale:
        the_radio.text = 'Don\'t Autoscale'
    else:
        the_radio.text = 'Autoscale'
    Start_Sim()

as_radio = button(text='Autoscale', pos=scene.title_anchor, bind=Auto_Scale)

def Allow_Unstable(the_radio):
    global allow_unstable, stop
    stop = True
    allow_unstable = not allow_unstable
    if allow_unstable:
        the_radio.text = 'Perfectly Aligned'
    else:
        the_radio.text = 'Imperfectly Aligned'
    Start_Sim()

au_radio = button(text='Imperfectly Aligned', pos=scene.title_anchor, bind=Allow_Unstable)

Start_Sim()    #start from default
