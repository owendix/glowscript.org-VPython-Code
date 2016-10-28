GlowScript 2.1 VPython
#This finds the location of the center of an underwater object at 
#different angles. It's always closer to the light entry point, in 
#both x and y, than the object.
#Here I can actually adjust the size, and use at least three converging light rays
#to find the location of two spots on the fish (front and back of ellipsoid, then
#scale for the cone-tail)


fish_color = color.cyan
water_color = color.cyan
water_opac = 0.25
img_color = color.yellow
img_opac = 0.75
ray_color = [color.yellow, color.red]
vray_opac = 0.5      #only for virtual ray
head_color=color.magenta
shirt_color=color.green
pants_color=color.orange

leftx = -2.5
#bigy = 2

canvw = 1000
canvh = (0.7)*canvw
scn = canvas(width=canvw, height=canvh, forward=vec(0,-0.125,-1),
            center=vec(0.39*leftx,0,0),range=-0.42*leftx, autoscale=False, 
            userzoom=False, userspin=False,up=vec(0,1,0)) #these are good vals

n_air = 1.
n_water = 4/3
n_wa = n_water/n_air

#Fish depth
d = -1 #meters
obsx0 = 1.25*d #not less (more negative) than than -2.25, not >= 0
Navg = 1        #1 gives accurate results to within what the human eye can detect

#Draw water level: pos = middle of box
waterlevel = box(pos=vec(0.5*leftx,0.6*d,0),
        length=3*abs(leftx),height=-1.2*d,width=3,color=water_color,
        opacity=water_opac)
    
#Draw fish object, light blue
f_L = 0.3       #half the fishes length
f_hw = 0.333*f_L
f_ax = vec(-f_L,0,0)
f_hw2 = f_hw/2
f_dtail = -0.666*f_ax
ray_r = 0.025*f_hw
fish = [ellipsoid(pos=vec(f_L/2,d,0),axis=f_ax,height=f_hw,width=f_hw,
        color=fish_color,up=vec(0,0,f_hw))]      #height points along up, out of screen
fish[0].nose = fish[0].pos + 0.5*fish[0].axis
fish[0].rear = fish[0].pos - 0.5*fish[0].axis

fish.append(cone(pos=fish[0].pos+f_dtail, axis=vec(-f_hw,0,0), 
                radius=f_hw2, color=fish_color))
        #[body,tail]
fish[1].fin = fish[1].pos + vec(0,fish[1].radius,0) #top of fin, for radius
#still shift by f_dtail (proportioned to length) along axis, now use 
#image's fin position to fine radius and height of the fish

instructs = label(text='Drag person around to see fish\'s image location',
            pos=fish[0].nose + vec(-2*f_L,0,0), height=18, box=False, line=False,
            opacity=0, align='right',font='monospace',color=color.red)

#make cylinders for rays
#ray1=sphere(pos=fish.pos,radius=0.1*fish.width,color=color.white,
#    make_trail=True,trail_radius=0.1*fish.width,opacity=0.2)#retain=10)
#ray2=sphere(pos=fish.pos,radius=0.1*fish.width,color=color.green,
#    make_trail=True,trail_radius=0.1*fish.width,opacity=0.2)#,retain=10)

#critical angle in radians
ang_crit = asin(1/n_wa)
#dang = (ang_crit-2*ang_w)/5

#scale dang for color scheme
#dcolor = dang/(ang_crit-ang_w)

#make observer
obs_r = f_hw2        #radius of person's head
obs_h = 4*f_hw   #height of person to eyes
obs_hr = obs_h-obs_r

#head first then body
#start observer's eye at x = 0
#light rays enter eye at front of head: observer[0].pos+observer[0].radius
observer = [sphere(pos=vec(obsx0-obs_r,obs_h,0),radius=obs_r,color=head_color)]
observer[0].eye = observer[0].pos+vec(observer[0].radius,0,0)
#when observer gets moved, move all parts including the eye
observer.append(cylinder(pos=observer[0].pos - vec(0,obs_r + 0.5*obs_hr,0),
                axis=vec(0,0.5*obs_hr,0),radius=obs_r,color=shirt_color))
observer.append(cylinder(pos=observer[1].pos - vec(0,0.5*obs_hr,0),
                axis=vec(0,0.5*obs_hr,0),radius=obs_r,color=pants_color))


#need to check angles for critical angle on-screen
#for optimization problem, solving other angle
def f(x,fish_pt):
    #assumes water level at y = 0
    #returns horizontal distance to left of fish_pt, where light will 
    #x is positive if to left, negative if to right
    #hit the water-boundary (y=0) as it moves toward the observer's eye
    
    sang_water = x/sqrt(x**2 + fish_pt.y**2)    #sin(ang_w)
    obs = observer[0].eye
    dx = (fish_pt.x - obs.x) - x    #positive or negative
    sang_air = dx/sqrt(dx**2 + obs.y**2)
    #return obs_x + obs_h*sin(x)/sqrt(n_wa**-2 - sin(x)**2) - f_L/2 - d*tan(x)
    return n_water*sang_water - n_air*sang_air
    
def bisect(func, param, low, high):
    #Find root of continuous function where f(low) and f(high) have opp signs
    #condition: low < high, and f(low) opp sign of f(high)
    #call with bisect(f, param, -f_L, fish_pt.x-observer[0].eye.x), always works
    tol = d*0.5e-5  #gives absolute scale for problem
    Nmax = 50
    #assert not samesign(func(low), func(high))

    for i in range(Nmax):
        midpoint = (low + high) / 2.0
        fm = func(midpoint, param)
        if fm == 0 or (high-low)/2 < tol:
            return midpoint
        if func(low, param)*fm> 0:     #same sign
            low = midpoint
        else:
            high = midpoint

    return midpoint

image=[]
w_rays=[]
a_rays=[]
v_rays=[]

#need to make this repeatedly callable
def make_image(is_initial):
    global image, w_rays, a_rays, v_rays
    #where light ray hits water: negative
    
    #take observer location and fish point location to find angle with bisect
    impact_x = [0,0,0]
    img = [0,0,0]
    #best way: invisible one to position of tail
    #visible one at top of fin
    for i, fish_pt in enumerate([fish[0].nose, fish[1].fin, fish[0].rear]):
        diffx = fish_pt.x - observer[0].eye.x
        xd = bisect(f,fish_pt,0,diffx)
        
        ang_w = atan(-xd/fish_pt.y)     #positive xd when d negative
                                #ang always between 0 and 90
        #xd = -d*tan(ang_w)    #positive, d neg, fish at (0,d,0)
        impact_x[i] = fish_pt.x-xd          #other rays shifted by s
        
        #perform shift(s)
        s = 0.001       #shift by 1 millimeter each time, find intersect --> img
        img[i] = vec(0,0,0)
        for n in range(Navg):
            s = (n+1)*s
            ang_wshift = atan((s+xd)/-fish_pt.y)
            ang_wimg = asin(n_wa*sin(ang_w))        #equals ang_air image
            ang_shift_wimg = asin(n_wa*sin(ang_wshift))     #equals angle_shift air image
            #calculate image location
            m = -1/tan(ang_wimg)
            m_s = -1/tan(ang_shift_wimg)  
            unshiftx = m_s*s/(m - m_s)      #positive distance from impact point
            img[i] = img[i] + vec(impact_x[i]+unshiftx, m*unshiftx, 0)
        img[i] = img[i]/Navg
    
    img_ax = img[0] - img[2]
    img_tailax = (f_hw/f_L)*img_ax
    img_h = mag(img_tailax)         #height points along up out of screen
    img_pos = img[0] - 0.5*img_ax
    img_dtail = proj(img[1]-img_pos, img_ax)        #assume fish is a straight line
    img_w = 2*mag(img_pos + img_dtail - img[1])
    if is_initial:
        image.append(ellipsoid(pos=img_pos,axis=img_ax,height=img_h,width=img_w,
                    up=vec(0,0,img_h),color=img_color,opacity=img_opac))
        image.append(cone(pos=image[0].pos+img_dtail, axis=img_tailax, 
                    radius=img_w/2, color=img_color,opacity=img_opac))
    
        #make light rays to the nose
        for n in range(2):
            impact = vec(impact_x[0] - n*s,0,0)
            w_rays.append(cylinder(pos=impact, radius=ray_r, 
                    axis=fish[0].nose-impact,color=ray_color[n]))
            a_rays.append(cylinder(pos=impact, radius=ray_r,
                    axis=observer[0].eye - impact,color=ray_color[n]))
            v_rays.append(cylinder(pos=impact, radius=ray_r,
                    axis=img[0]-impact,color=ray_color[n],
                    opacity=vray_opac))
        #make light rays ot the rear torso
        for n in range(2):
            impact = vec(impact_x[1] - n*s,0,0)
            w_rays.append(cylinder(pos=impact, radius=ray_r, 
                    axis=fish[1].fin-impact,color=ray_color[n]))
            a_rays.append(cylinder(pos=impact, radius=ray_r,
                    axis=observer[0].eye - impact,color=ray_color[n]))
            v_rays.append(cylinder(pos=impact, radius=ray_r,
                    axis=img[1]-impact,color=ray_color[n],
                    opacity=vray_opac))
    else:
        image[0].pos=img_pos
        image[0].axis=img_ax
        image[0].height=img_h
        image[0].width=img_w
        image[1].pos=image[0].pos+img_dtail
        image[1].axis=img_tailax
        image[1].radius=img_w/2
        
        for n in range(2):
            impact = vec(impact_x[0] - n*s,0,0)
            w_rays[n].pos=impact
            w_rays[n].axis=fish[0].nose-impact
            a_rays[n].pos=impact
            a_rays[n].axis=observer[0].eye-impact
            v_rays[n].pos=impact
            v_rays[n].axis=img[0]-impact
        for n in range(2,4):
            impact = vec(impact_x[1] - n*s,0,0)
            w_rays[n].pos=impact
            w_rays[n].axis=fish[1].fin-impact
            a_rays[n].pos=impact
            a_rays[n].axis=observer[0].eye-impact
            v_rays[n].pos=impact
            v_rays[n].axis=img[1]-impact

           
make_image(True)     #is_initial, don't do anything with this feature yet


#machinery for moving person (observer), cannot move fish
drag = False

def inside_object(a_pos):
    #approximate person as a rectangle, vertically aligned
    #position components of [right, top, left, bottom]
    wd2 = observer[0].radius
    rect=[observer[0].pos.x + wd2, 
        observer[0].pos.y + wd2,
        observer[0].pos.x - wd2,
        observer[-1].pos.y]
    if a_pos.x >= rect[0] or a_pos.x <= rect[2] or \
            a_pos.y >= rect[1] or a_pos.y <= rect[3]:
        return False
        
    return True

diff_vec0 = vec(0,0,0)

scn.bind("mousedown", def ():
    global drag, diff_vec0
    drag = inside_object(scn.mouse.pos)
    diff_vec0 = observer[0].pos - scn.mouse.pos   
)

scn.bind("mousemove", def ():
    global drag, observer
    eps=0.1*observer[0].radius
    if drag and observer[0].eye.x<-eps and observer[0].eye.y>observer[0].radius:
        new_pos = scn.mouse.pos + diff_vec0
        #check the y condition too
        if new_pos.x < -(observer[0].radius + eps) and new_pos.y > observer[0].radius:
            do1 = observer[1].pos - observer[0].pos
            do2 = observer[2].pos - observer[0].pos
            
            observer[0].pos = new_pos
            observer[0].eye = new_pos + vec(observer[0].radius,0,0)
            observer[1].pos = new_pos + do1
            observer[2].pos = new_pos + do2
        #if observer[0].pos.x > -(observer[0].radius+eps):
        #    for o in observer:
        #        o.pos = o.pos - vec(eps,0,0)
        #    observer[0].eye = observer[0].eye - vec(eps,0,0)
        #    drag = False
        
        make_image(False)
)
scn.bind("mouseup", def ():
    nonlocal drag
    drag = False
)
