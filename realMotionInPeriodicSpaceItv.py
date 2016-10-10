GlowScript 2.1 VPython
#spaceship with no retro rockets
#bullets assumed to be negligible mass compared to ship

the_rate = 1000
dt = 2**-8
Lfac = 15           #half of how much bigger scene is than ship
t = 0

Nworlds = 3
if Nworlds > 5:
    Nworlds = 5
if Nworlds < 0:
    Nworlds = 0

r_L = 37.2         #m, length of space shuttle endeavor
r_r = 23.8/2        #m, half wingspan of space shuttle endeavor
#the rest of this motion is not to scale with the endeavor
L = Lfac*r_L
world_min_r = r_L   #could also choose r_r, or (r_L+r_r)/2

a = 10*r_L           #m/s, constant linear acceleration
alpha = 3*pi         #rad/s^2, constant angular acceleration

v_max = 3*r_L       #m/s, max factor of ship length traveled
w_max = pi/4        #radians/second   #maximum angular speed

stabilize_rotation = True
degfac = 3600      #21600 rounds to nearest arcminute, 3600 rounds to deg/10

Nbullets = 20       #maximum number of bullets visible at one time
ibullet = 0         #mod this over Nbullet
vbullet = v_max     #can be changed
tbullet = 1.55*L/v_max       #time to travel 1.5*L if going at v_max, don't change
tbullet_lab = 2*tbullet/Nbullets
taccel_last = 0

lose = False

ship_color = color.red
bullet_color = color.white
fuel_color=color.green
canv_w = 700
canv_h = canv_w


scn_range = L   #I may need a factor
#create the boundaries: -L/L
scene=canvas(width=canv_w,height=canv_h,dir=vec(0,0,-1),
        range=scn_range,autoscale=False,userzoom=False,userspin=False,
        up=vec(0,1,0))

print('To navigate: left, up, right arrows')
print('To shoot: spacebar')

#function for periodic motion
def check_boundaries():
    global ship, bullets
    neg_diff = ship.pos - (-L*vec(1,1,0))
    pos_diff = ship.pos - (L*vec(1,1,0))
    if neg_diff.x < 0:
        ship.pos.x = L + neg_diff.x
    elif pos_diff.x > 0:
        ship.pos.x = -L + pos_diff.x
    if neg_diff.y < 0:
        ship.pos.y = L + neg_diff.y
    elif pos_diff.y > 0:
        ship.pos.y = -L + pos_diff.y
    for b in bullets:
        if b.visible:
            neg_diff = b.pos - (-L*vec(1,1,0))
            pos_diff = b.pos - (L*vec(1,1,0))
            if neg_diff.x < 0:
                b.pos.x = L + neg_diff.x
            elif pos_diff.x > 0:
                b.pos.x = -L + pos_diff.x
            if neg_diff.y < 0:
                b.pos.y = L + neg_diff.y
            elif pos_diff.y > 0:
                b.pos.y = -L + pos_diff.y
    
worlds = []
def make_worlds():
    wr_fac = 0.0625*L*(1/(1+Nworlds)+1)
    for n in range(Nworlds):
        rvec = vec.random()
        world_r = abs(rvec.z+2)*wr_fac
        if world_r < world_min_r:
            world_r = world_min_r   #must be set below max for good range
        rvec.z = 0
        rvec = rvec*(L-world_r)
        if mag(rvec) < world_r:
            rvec = rvec + (2*r_L+world_r)*hat(rvec) #don't put it in the center
        cvec = 0.5*(vec.random() + vec(1,1,1))
        for w in worlds:
            if mag(rvec - w.pos) < (world_r + w.radius):
                cvec = w.color      #use same color if they overlap
        worlds.append(sphere(pos=rvec,radius=world_r, 
            color=cvec))

#create the rocket (a cone)
ship = cone(pos=vec(0,0,0),axis=vec(0,r_L,0),radius=r_r,color=ship_color)
ship.v = vec(0,0,0)         #only along x or y
ship.w = vec(0,0,0)     #only along z or -z

bullets = []
def make_bullets():
    b_r = r_r/4
    for n in range(Nbullets):
        bullets.append(sphere(pos=vec(0,0,0), radius=b_r,
            color=bullet_color,visible=False))
        bullets[-1].v = vec(0,0,0)  #no velocity until they're fired
        bullets[-1].t = 0           #time until disappears: needed more than 
            #distance so it can't linger forever. Needs to be reused
        #I won't be updating their position until fired anyway

make_bullets()

#centroid of ship: solid cone its 1/4 the distance from base to apex r_L/4
#make it spin around centroid, may need rate(the_rate, move)
#want it to constantly be moving
if Nworlds > 0:
    make_worlds()
    
def check_collisions():
    #simple collision, if within the ships radius from sphere surface
    global lose
    for w in worlds:
        if w.visible:
            if mag(ship.pos + 0.25*ship.axis - w.pos) < (1.1*ship.radius + w.radius):
                label(text='Game Over\nReload Page To Restart', pos=vec(0,0,0),
                    height=36, color=color.red, align='center')
                lose = True
            
def check_target_shot():
    global worlds, bullets, lose
    for b in bullets:
        if b.visible:
            #check if I got shot
            if mag(ship.pos + 0.25*ship.axis - b.pos) < (1.1*ship.radius + b.radius):
                label(text='Game Over\nReload Page To Restart', pos=vec(0,0,0),
                    height=36, color=color.red, align='center')
                lose = True
            for w in worlds:
                if w.visible:
                    #check if worlds got shot
                    if mag(b.pos - w.pos) <= (w.radius + b.radius):
                        b.visible=False
                        w.radius = 0.75*w.radius
                        if w.radius < world_min_r:
                            w.visible = False

last_key = ' '

def key_input(evt):
    global ship, bullets, ibullet, taccel_last
    s = evt.key
    if s == 'left':     #positive z angular acceleration
        if (t - taccel_last) >= dt:
            taccel_last = t     #don't increment speed faster than dt
            if ship.w.z <= w_max:
                ship.w = ship.w + vec(0,0,alpha*dt)
            
    elif s == 'right':  #negative z angular acceleration
        if (t - taccel_last) >= dt:
            taccel_last = t     #don't increment speed faster than dt
            if ship.w.z >= -w_max:
                ship.w = ship.w + vec(0,0,-alpha*dt)
            
    elif s == 'up':        #set linear acceleration along axis
        if (t - taccel_last) >= dt:
            taccel_last = t     #don't increment speed faster than dt
            ship_speed = mag(ship.v)
            if mag(ship.v) <= v_max or dot(ship.v,ship.axis) < 0:        
                ship.v = ship.v + a*hat(ship.axis)*dt
                if mag(ship.v) > v_max:
                    ship.v = v_max*hat(ship.v)  #this should do it
                
    elif s == ' ':  #holding this down will let me know how fast keys repeat
        #fire a bullet, on release
        last_b = bullets[(ibullet+Nbullets-1)%Nbullets]
        shoot = False
        if last_b.visible:
            if last_b.t >= tbullet_lab:
                shoot = True
        else:
            shoot = True
        if shoot:
            bullets[ibullet].pos=ship.pos + ship.axis
            #make it relative to ship's velocity, here's the cool part
            bullets[ibullet].v = ship.v + vbullet*hat(ship.axis)
            bullets[ibullet].t = 0  #set lifetime
            bullets[ibullet].visible = True
            ibullet = (ibullet+1)%Nbullets
        
        
scene.bind('keydown', key_input)

scene.caption='Vary the bullet speed\n'
def set_vbullet(r):
    global vbullet
    vbullet = r.value

slider(min=0.5*v_max,value=v_max,max=1.5*v_max,length=canv_w,bind=set_vbullet)

n = 0
while not lose: #go forever
    rate(the_rate,wait)
    ship.pos = ship.pos + ship.v*dt
    #rotate around centroid=center of mass for a solid cone
    if stabilize_rotation:
        ang = round(ship.w.z*dt*degfac)/degfac
    else:
        ang = ship.w.z*dt
    ship.rotate(angle=ang, axis=vec(0,0,1),origin=ship.pos+0.25*ship.axis)
    #I am rounding to the nearest degree tenth of a degree 
    #perhaps improving control
    for b in bullets:
        if b.visible:
            b.pos = b.pos + b.v*dt
            b.t += dt
            if b.t > tbullet:
                b.visible=False
    check_target_shot()
        
    check_collisions()
    check_boundaries()
    n += 1
    t = n*dt
