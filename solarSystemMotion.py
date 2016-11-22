GlowScript 2.2 VPython
#Simulation of solar system motion using Newtonian gravity and an explicit midpoint-type method
#algorithm. Allows you to view real satellite motion as it orbits a moving planet.
cw = 800
ch = cw
the_rate0 = 1000
the_rate = the_rate0
#2**-9 is not accurate enough for moons with forward euler
#2**-8 is a little choppy with explicit midpoint-type method
dt0 = 2**-8         
dt = dt0
dtshrinkfac = 2**0
rangefac = 1.2
fwdfac = 1.125
font_h = 18
font_yoffset = 0.1*font_h
#I got it to get the right zoom but it still looks horrible
#radii of the worlds in meters, including the sun (0) to neptune
#Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune
Rs = 695700         #km
Re = 6371.008       #km
Rm = 1737.4         #km

lbl = ['Sun','Mercury','Venus',
        ['Earth','Moon'],
        ['Mars','Phobos','Deimos'],
        ['Jupiter','Io', 'Europa','Ganymede','Callisto'],
        ['Saturn','Mimas','Enceladus','Tethys','Dione','Rhea','Titan','Iapetus'],
        ['Uranus','Miranda','Ariel','Umbriel','Titania','Oberon'],
        ['Neptune', 'Triton']
        ]
clr = [vec(1,0.9,0.5), vec(0.9,0.9,1), vec(0.9,0.8,0.4),
        [vec(0.5,1,1), vec(0.9,0.9,1)],
        [vec(233/255,116/255,81/255),vec(233/255,175/255,125/255),vec(233/255,200/255,150/255)],
        [vec(233/255,150/255,100/255), vec(1,0.9,0.5), vec(233/255,200/255,150/255),color.gray(0.333), vec(0.5,0.25,0.25)],
        [vec(1,0.95,0.75),color.gray(0.25), color.gray(0.9), color.gray(0.8), color.gray(0.85), vec(1,0.975,0.85), vec(0,0.5,0.5),vec(1,0.975,0.9)],
        [vec(0,0.5,1), color.gray(0.85), vec(233/255,200/255,150/255), color.gray(0.5), color.gray(0.9), vec(1,0.9,0.9)],
        [vec(0.7,1,1), vec(1,0.8,0.8)]
        ]

#masses
#ms
sec2day = 86400
Ms = 1988500e24     #kg
GMs = 132712e6      #km^3/s^2 
GMs *= sec2day**2    #km^3/day^2
me = 5.9724e24      #kg
mm = 0.07346e24     #kg
#masses in kg
m = [Ms,0.33011e24,4.8675e24,
    [me,mm],
    [0.642e24,10.6e15,2.4e15],
    [1898e24, 893.2e20, 480.0e20, 1481.9e20, 1075.9e20],
    [568e24, 0.379e20, 1.08e20, 6.18e20, 11.0e20, 23.1e20, 1345.5e20, 18.1e20],
    [86.8e24, 0.66e20, 12.9e20, 12.2e20, 34.2e20, 28.8e20],
    [102e24, 214e20]
    ]

#bulk radius in km
R = [Rs, 2439.7, 6051.8,
    [Re,Rm],
    [3389.5, 10.25, 5.55],
    [69911, 1821.5, 1560.8, 2631.2, 2410.3],
    [58232, 198.9, 252.1, 531.1, 561.9, 763.3, 2575.6, 734.6],
    [25362, 236, 579, 584, 789, 761],
    [24622, 1353.4]
    ]
    
#saturns rings
#glowscript has a ring, but its a torus, I need a flat disk with a hole

G = 6.67408e-11     #m^3/(kg s^2)
G /= 1e9
G *= sec2day**2     #km^3/(kg day^2)

ogn = vec(0,0,0)  #treat sun as fixed
vogn = vec(0,0,0)
#perihelion in kilometers
#when eccentricities are small and 
#lunar characteristics sparse, average radius and speed used
perih = [0, 46.00e6, 107.48e6,
        [147.09e6,0.3633e6],
        [206.6e6,9378,23459],
        [740.5e6, 421.8e3, 671.1e3, 1070.4e3, 1882.7e3],
        [1352.6e6, 185.52e3, 238.02e3, 294.66e3, 377.40e3, 527.04e3, 1221.83e3, 3561.3e3],
        [2741.3e6, 129900.0, 199900.0, 266000.0, 436300.0, 583500.0],
        [4444.5e6, 354.76e3]
        ]
#vmax (at perihelion) in km/s
#some lunar speeds are approximated assuming a circle
vperih = [0,58.98,35.26,
        [30.29,1.076],
        [26.50,2.138,1.351],
        [13.72,17.338, 13.743, 10.88, 8.204],
        [10.18, 14.316, 12.633, 11.351, 10.028, 8.484, 5.572, 3.265],
        [7.11, 6.68,5.77, 4.67, 3.64, 3.15],
        [5.50, 4.39]
        ]
#orbital tilt in degrees, rel to sun
orbtilt = [0,7.00,3.39,
            [0,5.145],
            [1.9,1.08,1.79],
            [1.3, 0.04, 0.47, 0.18, 0.19],
            [2.5, 1.53, 0.00, 1.86, 0.02, 0.35, 0.33, 14.72],
            [0.8, 4.34, 0.04, 0.13, 0.08, 0.07],
            [1.8, 157.345]
            ]
r = []
v = []
for i in range(len(R)):
    if len(R[i]) > 0:
        r.append([])
        v.append([])
        for j in range(len(R[i])):
            vperih[i][j] *= sec2day
            orbtilt[i][j] = radians(orbtilt[i][j])
            r[i].append(perih[i][j]*vec(cos(orbtilt[i][j]),0,sin(orbtilt[i][j])))
            v[i].append(vperih[i][j]*vec(0,1,0))
            
    else:
        vperih[i] *= sec2day
        #assume planets lie perfectly in ecliptic: radians(0) = 0
        if i==0:
            r.append(ogn)
            v.append(vogn)
        else:
            orbtilt[i] = radians(orbtilt[i])
            r.append(perih[i]*vec(cos(orbtilt[i]),0,sin(orbtilt[i])))
            v.append(vperih[i]*vec(0,1,0))

#Earth orbit
Te = 365.256            #Earth sidereal orbit in days, (sidereal: wrt fixed stars)

rlast = r[-1]
rmax = 0
if len(rlast) > 0:
    for rl in rlast:
        rmax += mag(rl)
else:
    rmax = mag(rlast)
scn = canvas(width=cw, height=ch, center=ogn, forward=vec(0,0,-1),up=vec(0,1,0),
            range=rangefac*rmax)
scn.center0 = scn.center
scn.range0 = scn.range
scn.forward0 = scn.forward
scn.camerapos0 = scn.camera.pos
scn.up0=scn.up

instr = label(pos=vec(0,-0.9*rmax,0),
            text='0-8: Focus on Sun/planets, l: Toggle labels, t: Toggle trails (improves speed), \n'
            +'Scroll: Zoom, Ctrl-drag: Reorient\n'+
            +'Uses Newtonian gravity and an explicit midpoint-type method\n'
            +'Some satellite motion approximated as circular', height=font_h, color=color.white,
            align='center', font='serif', box=False, line=False, opacity=0, visible=True)

#start the moon assuming it is at its closest point out and high in the 
#orbital tilt, traveling with Earth, so magnitudes add, and axial tilt is toward
#the sun, and the Earth is at its perihelion
vis_label = True

body = []
oldr = []
tmpr = []
oldv = []
tmpv = []
for i, R_i in enumerate(R):
    if len(R_i) > 0:
        r0 = r[0]
        v0 = v[0]
        bi = []
        oldr.append([])
        oldv.append([])
        tmpr.append([])
        tmpv.append([])
        for j, R_j in enumerate(R_i):
            bi.append(sphere(pos=r0+r[i][j], radius=R_j, color=clr[i][j], make_trail=True))
            oldr[i].append(bi[j].pos)
            tmpr[i].append(bi[j].pos)
            bi[j].lbl=label(text=lbl[i][j],pos=bi[j].pos+bi[j].radius*vec(0,1,0),
                    height=font_h,box=False,yoffset=font_yoffset*(j+1),
                    line=False,opacity=0,align='center',color=clr[i][j],
                    visible=(j==0 and vis_label))
            
            
            bi[j].v = v0+v[i][j]
            oldv[i].append(bi[j].v)
            tmpv[i].append(bi[j].v)
            bi[j].m = m[i][j]

            if j == 0:
                r0 = r[i][0]
                v0 = v[i][0]
        body.append(bi.copy())
    else:
        if i == 0:
            #need this because of addition, not ogn and vogn
            r0 = vec(0,0,0)
            v0 = vec(0,0,0)
        else:
            r0 = r[0]
            v0 = v[0]
        #print(r0,r[i],v0,R_i,clr[i])
        body.append(sphere(pos=r0+r[i], radius=R_i, color=clr[i], make_trail=(i!=0)))
        oldr.append(body[i].pos)
        tmpr.append(body[i].pos)
        body[i].lbl=label(text=lbl[i],pos=body[i].pos+body[i].radius*vec(0,1,0),
                    height=font_h,box=False,yoffset=font_yoffset,line=False,
                    opacity=0,color=clr[i],align='center',visible=vis_label)
        body[i].v = v0 + v[i]
        oldv.append(body[i].v)
        tmpv.append(body[i].pos)
        body[i].m = m[i]
        #if i==some number that is not the sun, add an axis, pov_loc, pov_fwd
        

def toggle_label():
    global body, vis_label
    
    #new status for if labels are visible
    vis_label = not body[0].lbl.visible
    
    for i, b in enumerate(body):
        if len(R[i]) > 0:
            for j in range(len(R[i])):
                if following[i]:
                    body[i][j].lbl.visible = vis_label
                else:
                    body[i][j].lbl.visible = (j==0 and vis_label)

        else:
            b.lbl.visible = vis_label
    

def toggle_trail():
    #toggle all trails on and off
    global body
        
    for i, b in enumerate(body):
        if len(R[i]) > 0:
            for B in b:
                B.make_trail = not B.make_trail
                if not B.make_trail:
                    B.clear_trail()
        else:
            if i != 0:
                b.make_trail = not b.make_trail
                if not b.make_trail:
                    b.clear_trail()
            

#calculations to change perspective
following = [False for x in range(len(R))]
#pov = [False for x in range(len(R))]
#pov_loc = sphere(pos=ogn,radius=1e-6,color=color.black,opacity=0)
#pov_fwd = sphere(pos=ogn,radius=1e-6,color=color.black,opacity=0)
shrink_dt = False
unshrink_dt = False
def change_view(evt):
    global following, pov, scn, shrink_dt, unshrink_dt, body
    #need to make it so pushing k turns everything else to False even if different k
    
    k = evt.key
    if k == 't':
        toggle_trail()
        return
    elif k =='l':
        toggle_label()
        return
    
    k = int(k)
    if k in range(len(R)):
        #print('k,following[k], pov[k]:',k,following[k],pov[k])
        
        if following[k]:#return to sun
            unshrink_dt = True
            scn.camera.follow(None)
            scn.center = scn.center0
            scn.range = scn.range0
            scn.forward = scn.forward0
            scn.camera.pos = scn.camerapos0
            scn.up = scn.up0
            
            if len(R[k]) > 0:
                #turn moon labels off
                for j in range(len(R[k])):
                    body[k][j].lbl.visible=(j==0 and vis_label)
            
            for i in range(len(R)):
                #pov[i] = False
                following[i] = False
            
        else:#go to following
            shrink_dt = k!=0
            if len(R[k]) > 0:
                scn.camera.follow(body[k][0])
            else:
                scn.camera.follow(body[k])
            scn.center = scn.center0
            scn.range = scn.range0
            scn.forward = scn.forward0
            scn.up = scn.up0
            
            for i in range(len(R)):
                if i == k:
                    following[i] = True
                    if len(R[i]) > 0:
                        #turn moon labels on
                        for j in range(len(R[i])):
                            body[i][j].lbl.visible=vis_label
                else:
                    following[i] = False
                    if len(R[i]) > 0:
                        #turn other moon labels off
                        for j in range(len(R[i])):
                            body[i][j].lbl.visible=(j==0 and vis_label)
                    
                
                #pov[i] = False
                

scn.bind('keyup',change_view)

#vary animation rate
def vary_rate(s):
    global the_rate
    the_rate = exp(s.value)

scn.caption='Vary the animation rate:\n'
slider(min=log(the_rate0/100), value=log(the_rate), max=log(the_rate0*100), length=cw, bind=vary_rate)

#iterate this in time
t = 0
n = 0
n_print = 100

while True:
    rate(the_rate)
    #adjust dt for midpoint method, first step
    if n == 0:
        dt /= 2
        #want positions for a to be current time
        #don't change current until the end, store updated in tmp
        #assign current to old and tmp to current
        #at n==0, oldr and oldv are body.pos and body.v
    

    for i, b in enumerate(body):
        if len(R[i]) > 0:
            for j, B in enumerate(b):
                #calculate the acceleration due to all other masses
                a = vec(0,0,0)
                for k, c in enumerate(body):
                    if len(R[k]) > 0:
                        for l, C in enumerate(c):
                            if k != i or l != j:
                                Gm = G*C.m
                                a = a + (Gm/mag2(C.pos - B.pos))*hat(C.pos - B.pos)
                    else:
                        if k == 0:
                            Gm = GMs
                        else:
                            Gm = G*c.m
                        a = a + (Gm/mag2(c.pos - B.pos))*hat(c.pos - B.pos)
                #use accel and update position of body
                tmpr[i][j] = oldr[i][j] + B.v*2*dt
                tmpv[i][j] = oldv[i][j] + a*2*dt
                B.lbl.pos=tmpr[i][j] + B.radius*vec(0,1,0)
                
        else:
            #don't update Sun's position: assume constant
            if i != 0:
                #calculate the acceleration due to all other masses
                a = vec(0,0,0)
                for k, c in enumerate(body):
                    if len(R[k]) > 0:
                        for l, C in enumerate(c):
                            Gm = G*C.m
                            a = a + (Gm/mag2(C.pos - b.pos))*hat(C.pos - b.pos)
                    else:
                        if k != i:
                            if k == 0:
                                Gm = GMs
                            else:
                                Gm = G*c.m
                            a = a + (Gm/mag2(c.pos - b.pos))*hat(c.pos - b.pos)
                
                #use accel and update position of body
                tmpr[i] = oldr[i] + b.v*2*dt
                tmpv[i] = oldv[i] + a*2*dt
                b.lbl.pos=tmpr[i] + b.radius*vec(0,1,0)
    
    #unadjust dt for midpoint method, the rest
    if n == 0:
        dt *= 2
    #copy tmp (current) to old then tmp to body
    for i in range(len(R)):
        if len(R[i]) > 0:
            for j in range(len(R[i])):
                oldr[i][j] = body[i][j].pos
                oldv[i][j] = body[i][j].v
                body[i][j].pos = tmpr[i][j]
                body[i][j].v = tmpv[i][j]
        else:
            oldr[i] = body[i].pos
            oldv[i] = body[i].v
            body[i].pos = tmpr[i]
            body[i].v = tmpv[i]
                
        
    n += 1
    t += dt
    
    if shrink_dt:
        dt = dt0*dtshrinkfac
        shrink_dt = False
    if unshrink_dt:
        dt = dt0
        unshrink_dt = False
    
    if t >= Te - 0.5*dt and t <= Te + 0.5*dt:
        for i, b in enumerate(body):
            if len(R[i]) > 0:
                for B in b:
                    B.clear_trail()
            else:
                b.clear_trail()
