GlowScript 2.1 VPython
#Drop three perfectly elastic balls of different sizes

auto_scale = False   #follows ball
allow_unstable = True

g = 9.8          #m/s^2, =0 and =9.8 to see effects (minimal)
dt = 1e-4
the_rate = 1e4
eps = 1e-10
#t_done = 20    #Currently running indefinitely

ms = [1, 0.5, 0.1]
rs = [0.1, 0.07937, 0.04642]
space = rs[0]/3
colors = [color.yellow, color.cyan, color.magenta] 
y_pos0 = 1.0        #m

canv_width=1000
canv_height=canv_width*(9/16)
scene=display(width=canv_width,height=canv_height,center=vec(0,y_pos0,0))
g_thick = rs[0]
scene.range=20*g_thick
scene.autoscale=auto_scale
g_zlim = 10*g_thick
g_xlim = 40*g_thick

ground = box(pos=vec(0,-0.5*g_thick,0),axis=vec(0,1,0),length=g_thick,
            width=2*g_zlim,height=2*g_xlim,color=color.green)

balls = []
xz_comps = vec(0,0,0)
for i in range(len(rs)):
    if i == 0:
        the_y = y_pos0
    else:
        the_y = old_b_pos.y + rs[i-1] + rs[i] + space
    if allow_unstable:
        xz_comps = (rs[2]*eps)*vec.random() #tiny random perturbation
    balls.append(sphere(pos=vec(xz_comps.x,the_y,xz_comps.z),radius=rs[i], 
        color=colors[i]))
    balls[-1].m = ms[i]
    balls[-1].v = vec(0,0,0)
    balls[-1].a = vec(0,-g,0)
    old_b_pos = balls[-1].pos

#no front wall

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
                if a.pos.y > b.pos.y:
                    print('t swap:',t,'i,j:',i,j+1,'a.m, b.m',a.m, b.m)

    return balls

#P = 0
t = 0
n = 0
dp = 0
y_max = balls[2].pos.y
while y_max > -g_thick:
    rate(the_rate)
    
    #make them move
    for b in balls:
        b.pos = b.pos + b.v*dt
        b.v = b.v + b.a*dt

    #check for ground collisions
    [balls,dp] = check_ground(balls,g_xlim,g_zlim)
    
    #check for interparticle collisions
    balls = collide_balls(balls)
    #if t == 0.53771:
    #    scene.pause()

    y_max = 0
    for b in balls:
        if b.pos.y > y_max:
            y_max = b.pos.y
    n += 1
    t = n*dt
