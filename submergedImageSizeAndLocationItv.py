GlowScript 2.1 VPython
#This finds the location of the center of an underwater object at 
#different angles. It's always closer to the light entry point, in 
#both x and y, than the object.
#Here I can actually adjust the size, and use at least three converging light rays
#to find the location of two spots on the fish (front and back of ellipsoid, then
#scale for the cone-tail)
#tail width, as a cone, is false

bg_color = color.gray(0.75)
fish_color = color.yellow
fishfin_color = color.red
water_color = color.blue
water_opac = 0.2
img_color = vec(0,1,0.25)        #fish_color #color.yellow
imgfin_color = color.magenta
img_opac = 0.8
ray_color = [color.yellow, color.red]
vray_opac = 0.5      #only for virtual ray
head_color=color.magenta
shirt_color=color.yellow
pants_color=color.cyan
silver = vec(0.75,0.75,0.75)
bubble_color = silver + vec(-0.125, 0, 0.125)
img_bubble_color = silver + vec(-0.25,0,0.25)

the_rate0 = 100     #used for moving a bubble
the_rate = the_rate0
dt = 2**-8
n_try_bubble = 100
p_bubble = 0.333

leftx = -2.5
#bigy = 2

canvw = 1100
canvh = (0.7)*canvw
the_range = -0.44*leftx
the_center=vec(0.39*leftx,0,0)
#the_forward=vec(0,-0.125,-1)
the_forward = vec(0,0,-1)
the_up=vec(0,1,0)       #I don't think I need to change this
scn = canvas(width=canvw, height=canvh, forward=the_forward,background=bg_color,
            center=the_center,range=the_range, autoscale=False, 
            userzoom=False, userspin=False,up=the_up) #these are good vals
scn.lights=[distant_light(direction=vec(0, 1, 0), color=color.gray(0.8)),
   distant_light(direction=vec(-1, 1, 0), color=color.gray(0.3))]
scn.ambient = color.gray(0.45)
#direction gives from which direction, not in which direction does it point
            
the_camerapos = scn.camera.pos      #get the 

scn.caption = "Vary the animation rate: \n"
def set_rate(s):
    global the_rate
    the_rate = s.value

slider(min=the_rate0/10,value=the_rate0,max=the_rate0*10,length=canvw,bind=set_rate)

n_air = 1.
n_water = 4/3
n_wa = n_water/n_air

#Fish depth
d = -1 #meters
s = 0.001   #slight shift of light rays to find intersection
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
tail_hfac = 0.25      #height factor for tail, height along up --> z
fin_hfac = 0.15      #height factor for tail, height along up
fin_lwfac = 0.5      #length,width factor for fins, not for use on tailfin
#fish body
fish = [ellipsoid(pos=vec(f_L/2,d,0),axis=f_ax,height=f_hw,width=f_hw,
        color=fish_color,up=vec(0,0,f_hw))]      #height points along up, out of screen
fish[0].nose = fish[0].pos + 0.5*fish[0].axis
fish[0].rear = fish[0].pos - 0.5*fish[0].axis
fish[0].flank = fish[0].pos + 0.5*fish[0].up
#fish tail: 1
fish.append(pyramid(pos=fish[0].pos+f_dtail, axis=vec(-f_hw,0,0), 
                width=f_hw, height=tail_hfac*f_hw, up=vec(0,0,tail_hfac*f_hw),
                color=fishfin_color))
        #[body,tail]
#fish[1].fin = fish[1].pos
fish[1].fintip = fish[1].pos + vec(0,fish[1].width/2,0) #top of fin, for radius
#fish eyes: 2, 3
fish.append(sphere(pos=fish[0].pos + 0.4*f_ax+0.2*(fish[0].up+cross(hat(f_ax),fish[0].up)), 
                radius=0.1*fish[0].height, color=color.black))
fish.append(sphere(pos=fish[0].pos + 0.4*f_ax+0.2*(-fish[0].up+cross(hat(f_ax),fish[0].up)), 
                radius=0.1*fish[0].height, color=color.black))
fish[0].d_dorsal = f_hw2*hat(cross(fish[0].axis,fish[0].up))
fish[0].dorsal = fish[0].pos + fish[0].d_dorsal
#top dorsal fin: 4
fish.append(box(pos=fish[0].pos+0.95*fish[0].d_dorsal,
            axis=fin_lwfac*(2*fish[0].d_dorsal-f_hw2*hat(fish[0].axis)), 
            width=fin_lwfac*f_hw,height=fin_hfac*f_hw, 
            up=vec(0,0,fin_hfac*f_hw),color=fishfin_color))
#bottom dorsal fin: 5
fish.append(box(pos=fish[0].pos-0.95*fish[0].d_dorsal,
            axis=-fin_lwfac*(2*fish[0].d_dorsal+f_hw2*hat(fish[0].axis)), 
            width=fin_lwfac*f_hw,height=fin_hfac*f_hw, 
            up=vec(0,0,fin_hfac*f_hw),color=fishfin_color))
#left pectoral fin: 6
fish.append(pyramid(pos=fish[0].pos+0.8*fish[0].up+0.2*f_ax,
            axis=-fin_lwfac*fish[0].up, width=fin_lwfac*f_hw,
            height=fin_hfac*f_hw, up=fin_hfac*f_hw*hat(fish[0].axis),
            color=fishfin_color))
#right pectoral fin: 7
fish.append(pyramid(pos=fish[0].pos-0.8*fish[0].up+0.2*f_ax,
            axis=fin_lwfac*fish[0].up, width=fin_lwfac*f_hw,
            height=fin_hfac*f_hw, up=fin_hfac*f_hw*hat(fish[0].axis),
            color=fishfin_color))
                
#still shift by f_dtail (proportioned to length) along axis, now use 
#image's fin position to find radius and height of the fish

instructs = label(text='Drag the person around to see the fish\'s image size and location.\n'+
            'An image point is placed where virtual light rays converge.\n'+
            'Press spacebar to see from the person\'s eyes.\n'+
            'Fish\'s image is approximated as piecewise linear.\n'+
            '\tLake bottom edges are unadjusted. Bubble physics is ignored.\n'+
            'Bubble images are approximated as spheres. Colors are false.',
            pos=fish[0].nose + vec(-1.05*f_L,2.8*f_hw,0), height=18, box=False, line=False,
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
def f(x,obj_pt, is_z0):
    #assumes water level at y = 0
    #returns horizontal distance to left of fish_pt, where light will 
    #x is positive if to left, negative if to right
    #hit the water-boundary (y=0) as it moves toward the observer's eye
    
    sang_water = x/sqrt(x**2 + obj_pt.y**2)    #sin(ang_w)
    obs = observer[0].eye
    if not is_z0:
        obj_obs=obj_pt - obs
        obj_obs.y = 0
        dx = mag(obj_obs) - x               #can't pass fish observe point
    else:
        dx = (obj_pt.x - obs.x) - x    #positive or negative
    sang_air = dx/sqrt(dx**2 + obs.y**2)
    #return obs_x + obs_h*sin(x)/sqrt(n_wa**-2 - sin(x)**2) - f_L/2 - d*tan(x)
    return n_water*sang_water - n_air*sang_air
    
def bisect(func, param1, param2, low, high):
    #Find root of continuous function where f(low) and f(high) have opp signs
    #condition: low < high, and f(low) opp sign of f(high)
    #call with bisect(f, param, -f_L, fish_pt.x-observer[0].eye.x), always works
    tol = d*0.5e-5  #gives absolute scale for problem
    Nmax = 50
    #assert not samesign(func(low), func(high))

    for i in range(Nmax):
        midpoint = (low + high) / 2.0
        fm = func(midpoint, param1, param2)
        if fm == 0 or (high-low)/2 < tol:
            return midpoint
        if func(low, param1, param2)*fm> 0:     #same sign
            low = midpoint
        else:
            high = midpoint

    return midpoint
    
def get_image_points(obj_pt, is_z0):
    #convert to observer-based coordinate system in case object out of z=0 plane
    global s
    if not is_z0:
        obj_obs = obj_pt - observer[0].eye
        obj_obs.y = 0
        
        diffx = mag(obj_obs)
        #diffx = obj_pt.x - observer[0].eye.x
        xd = bisect(f,obj_pt,is_z0,0,diffx)
        impact_x = diffx - xd
        
        ang_w = atan(-xd/obj_pt.y)     #positive xd when d negative
                                #ang always between 0 and 90
        #xd = -d*tan(ang_w)    #positive, d neg, fish at (0,d,0)
        #impact_x[i] = obj_pt.x-xd          #other rays shifted by s
        
        #perform shift(s)
        s = 0.001       #shift by 1 millimeter each time, find intersect --> img
        img = vec(0,0,0)
        for n in range(Navg):
            s = (n+1)*s
            ang_wshift = atan((s+xd)/-obj_pt.y)
            ang_wimg = asin(n_wa*sin(ang_w))        #equals ang_air image
            ang_shift_wimg = asin(n_wa*sin(ang_wshift))     #equals angle_shift air image
            #calculate image location
            m = -1/tan(ang_wimg)
            m_s = -1/tan(ang_shift_wimg)  
            unshiftx = m_s*s/(m - m_s)      #positive distance from impact point
            img = img + vec(impact_x+unshiftx, m*unshiftx, 0)
        img = img/Navg
        #convert back from observer based coord system
        u = hat(obj_obs)
        img = vec(img.x*u.x, img.y, img.x*u.z)
        impact_x = vec(impact_x*u.x, 0, impact_x*u.z)
        obs_0 = vec(observer[0].eye.x,0,observer[0].eye.z)
        img = img + obs_0
        impact_x = impact_x + obs_0  #not needed w/o rays
    else:
        diffx = obj_pt.x - observer[0].eye.x
        xd = bisect(f,obj_pt,is_z0,0,diffx)
        
        ang_w = atan(-xd/obj_pt.y)     #positive xd when d negative
                                #ang always between 0 and 90
        #xd = -d*tan(ang_w)    #positive, d neg, fish at (0,d,0)
        impact_x = obj_pt.x-xd          #other rays shifted by s
        
        #perform shift(s)
        s = 0.001       #shift by 1 millimeter each time, find intersect --> img
        img = vec(0,0,0)
        for n in range(Navg):
            s = (n+1)*s
            ang_wshift = atan((s+xd)/-obj_pt.y)
            ang_wimg = asin(n_wa*sin(ang_w))        #equals ang_air image
            ang_shift_wimg = asin(n_wa*sin(ang_wshift))     #equals angle_shift air image
            #calculate image location
            m = -1/tan(ang_wimg)
            m_s = -1/tan(ang_shift_wimg)  
            unshiftx = m_s*s/(m - m_s)      #positive distance from impact point
            img = img + vec(impact_x+unshiftx, m*unshiftx, 0)
        img = img/Navg
        
    return [img, impact_x]
    

image=[]
w_rays=[]
a_rays=[]
v_rays=[]
bubbles = []
image_bubbles = []

#need to make this repeatedly callable
def make_image(is_initial):
    global image, w_rays, a_rays, v_rays
    #where light ray hits water: negative
    
    #take observer location and fish point location to find angle with bisect
    impact_x = [0,0,0,0,0]
    img = [0,0,0,0,0]
    #best way: invisible one to position of tail
    #visible one at top of fin
    for i, fish_pt in enumerate([fish[0].nose, fish[1].fintip, fish[0].rear, 
            fish[0].dorsal,fish[1].pos,fish[0].flank]):#4 is tail pos
        if i < 5:
            [img[i],impact_x[i]] = get_image_points(fish_pt,True)
        else:
            [img[i],impact_x[i]] = get_image_points(fish_pt,False)
        
    img_ax = img[0] - img[2]
    img_pos = img[0] - 0.5*img_ax
    img_dtail = proj(img[1]-img_pos, img_ax)    #assume fish body is a straight line
    #fish flank
    tmpvec = img[5] - img[0]
    #img_h =(f_hw/f_L)*mag(img_ax)         #height points along up out of screen
    img_h = 2*mag(tmpvec - proj(tmpvec,-img_ax))
    img_up = vec(0,0,img_h)
    
    #get img_w2 from dorsal: img[3], assume its directly in middle of torso
    tmpvec = img[3]-img[0]
    img_w2 = mag(tmpvec - proj(tmpvec,-img_ax))
    #img_w2 = mag(img_pos + img_dtail - img[1])
    img_w = 2*img_w2
    tmpvec = img[4]-img[1]      #from tail fintip img to tail fin pos img
    #vector pointing to shifted tail pos from fintip
    img_taild_pos = tmpvec - proj(tmpvec,img[2]-img[4])  
    img_tailpos = img[1] + img_taild_pos
    #half of tail pyramid outside body (factor of 2)
    img_tailax = 2*(img[2]-img_tailpos)
    img_tailw = 2*mag(img_taild_pos)
    img_tailh = tail_hfac*(f_hw/f_L)*mag(img_tailax)

    img_lpecpos = img_pos+0.8*img_up+0.2*img_ax
    img_rpecpos = img_pos-0.8*img_up+0.2*img_ax
    img_lpecax = -0.5*img_up
    img_pech = fin_hfac*img_h
    img_pecw = fin_lwfac*img_w
    img_pecup = fin_hfac*img_h*hat(img_ax)

    if is_initial:
        #image body: 0
        image.append(ellipsoid(pos=img_pos,axis=img_ax,height=img_h,width=img_w,
                    up=img_up,color=img_color,opacity=img_opac))
        #need image of fin tip to line up properly
        #image tail: 1
        image.append(pyramid(pos=img_tailpos, axis=img_tailax, 
                    width=img_tailw, height=img_tailh, up=vec(0,0,img_tailh),
                    color=imgfin_color,opacity=img_opac))
        #image eyes: 2, 3
        image.append(sphere(pos=image[0].pos + 0.4*img_ax+
                0.2*(image[0].up+img_w*cross(hat(img_ax),hat(image[0].up))), 
                radius=0.1*image[0].height, color=color.black))
        image.append(sphere(pos=image[0].pos + 0.4*img_ax+
                0.2*(-image[0].up+img_w*cross(hat(img_ax),hat(image[0].up))), 
                radius=0.1*image[0].height, color=color.black))
        image[0].d_dorsal = img_w2*hat(cross(image[0].axis,image[0].up))
        #top dorsal fin: 4
        image.append(box(pos=image[0].pos+0.95*image[0].d_dorsal,
                axis=fin_lwfac*(2*image[0].d_dorsal-img_w2*hat(image[0].axis)), 
                width=fin_lwfac*img_h,height=fin_hfac*img_h, 
                up=vec(0,0,fin_hfac*img_h),color=imgfin_color, opacity=img_opac))
        #bottom dorsal fin: 5
        image.append(box(pos=image[0].pos-0.95*image[0].d_dorsal,
                axis=-fin_lwfac*(2*image[0].d_dorsal+img_w2*hat(image[0].axis)), 
                width=fin_lwfac*img_h,height=fin_hfac*img_h, 
                up=vec(0,0,fin_hfac*img_h),color=imgfin_color, opacity=img_opac))
        #left pectoral fin: 6
        image.append(pyramid(pos=img_lpecpos,
                axis=img_lpecax,width=img_pecw, 
                height=img_pech, up=img_pecup,color=imgfin_color, 
                opacity=img_opac))
        #right pectoral fin: 7
        image.append(pyramid(pos=img_rpecpos,
                axis=-img_lpecax, width=img_pecw,
                height=img_pech, up=img_pecup,
                color=imgfin_color, opacity=img_opac))
    
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
                    axis=fish[1].fintip-impact,color=ray_color[n]))
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
        image[0].up = img_h*vec(0,0,1)
        image[1].pos=img_tailpos
        image[1].axis=img_tailax
        image[1].height=img_tailh
        image[1].width=img_tailw
        image[2].pos=image[0].pos + 0.4*img_ax+ \
                0.2*(image[0].up+img_w*hat(cross(img_ax,image[0].up)))
        image[2].radius=0.1*image[0].height
        image[3].pos=image[0].pos + 0.4*img_ax+ \
                0.2*(-image[0].up+img_w*hat(cross(img_ax,image[0].up)))
        image[3].radius=0.1*image[0].height
        #4,5,6,7 are dorsal, then lateral fins
        image[0].dorsal = img_w2*hat(cross(image[0].axis,image[0].up))
        image[4].pos=image[0].pos+0.95*image[0].dorsal
        image[4].axis=fin_lwfac*(2*image[0].dorsal-img_w2*hat(image[0].axis)) 
        image[4].width=fin_lwfac*img_h
        image[4].height=fin_hfac*img_h
        image[5].pos=image[0].pos-0.95*image[0].dorsal
        image[5].axis=-fin_lwfac*(2*image[0].dorsal+img_w2*hat(image[0].axis)) 
        image[5].width=fin_lwfac*img_h
        image[5].height=fin_hfac*img_h
        image[6].pos=img_lpecpos
        image[6].axis=img_lpecax
        image[6].width=img_pecw
        image[6].height=img_pech
        image[6].up=img_pecup
        image[7].pos=img_rpecpos
        image[7].axis=-img_lpecax
        image[7].width=img_pecw
        image[7].height=img_pech
        image[7].up=img_pecup
        
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
            w_rays[n].axis=fish[1].fintip-impact
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
    if not pov:
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
        
        make_image(False)
)

scn.bind("mouseup", def ():
    nonlocal drag
    drag = False
)

pov = False

def key_input(evt):
    global pov, a_rays, w_rays, v_rays
    s = evt.key
    if s == ' ':
        if not pov:
            scn.center = image[0].pos
            scn.forward = hat((a_rays[0].axis + a_rays[2].axis)/(-2))
            scn.camera.pos = observer[0].eye
            scn.up = cross(scn.forward,vec(0,0,-1))

        else:
            scn.center = the_center
            scn.forward = the_forward
            scn.camera.pos = the_camerapos
            scn.up = the_up

        pov = not pov
        for f in fish:
            f.visible = not f.visible
        instructs.visible = not instructs.visible
        for rays in [a_rays, w_rays, v_rays]:
            for r in rays:
                r.visible = not r.visible
    

#wait for key to be released to make move, so it doesn't switch really fast
scn.bind('keyup', key_input)

#every now and then make a bubble appear, along with its image, don't show rays
t = 0
n_t = 0
Tfin = 1
ang_vfac = 2    #max radians per second
bubble_pos0 = 0.5*(fish[0].nose + fish[0].pos) #+ fish[0].d_dorsal)#now gills bubble
b_v = vec(0,1,0)         #m/s
while True:
    rate(the_rate)
    
    
    ang_v1 = ang_vfac*sin(2*pi*t/Tfin)     #radians/sec
    ang_v2 = ang_vfac*sin(pi*(2*t/Tfin + 0.5))
    
    fish[6].rotate(angle=ang_v1*dt, axis=vec(0,1,0), 
                    origin=fish[6].pos + 0.625*fish[6].axis)
    fish[6].rotate(angle=ang_v2*dt, axis=-fish[6].axis,
                    origin=fish[6].pos + fish[6].axis)
    fish[7].rotate(angle=ang_v1*dt, axis=vec(0,-1,0),
                    origin=fish[7].pos + 0.625*fish[7].axis)
    fish[7].rotate(angle=ang_v2*dt, axis=fish[7].axis,
                    origin=fish[7].pos + fish[7].axis)
    dorsal = hat(cross(image[0].axis,image[0].up))
    image[6].rotate(angle=ang_v1*dt, axis=dorsal, 
                    origin=image[6].pos + 0.625*image[6].axis)
    image[6].rotate(angle=ang_v2*dt, axis=-image[6].axis,
                    origin=image[6].pos + image[6].axis)
    image[7].rotate(angle=ang_v1*dt, axis=-dorsal,
                    origin=image[7].pos + 0.625*image[7].axis)
    image[7].rotate(angle=ang_v2*dt, axis=image[7].axis,
                    origin=image[7].pos + image[7].axis)
    
    
    if n_t%n_try_bubble == 0:
        v = vec.random()
        if abs(v.x) > p_bubble:
            #create a bubble, at gill level ~halfway between pos and nose
            bbl_r = (0.5 + abs(v.y))*fish[2].radius     #0.5-1.5 times fish eye radius
            if v.z > 0:
                bubble_posi = bubble_pos0 + 0.5*fish[0].up      #shift for gills
            else:
                bubble_posi = bubble_pos0 - 0.5*fish[0].up
            #calculate image bubble position and size, image location for it 
            #and bottom of bubble
            b_impact_x = [0,0]
            b_img = [0,0]
            for i, obj_pt in enumerate([bubble_posi,bubble_posi+vec(0,bbl_r,0)]):
                [b_img[i],b_impact_x[i]] = get_image_points(obj_pt,False)
                
            #create bubble and image of bubble
            #image calculation assumes its in the z = 0 plane
            b_img_r = mag(b_img[1]-b_img[0])
            create_bubble = True
            for i,b in enumerate(bubbles):
                if b.popped:
                    b.pos=bubble_posi
                    b.radius=bbl_r
                    b.popped = False
                    b.visible=True
                    #image bubble
                    image_bubbles[i].pos=b_img[1] - vec(0,b_img_r,0)
                    image_bubbles[i].radius=b_img_r
                    image_bubbles[i].visible=True
                    create_bubble = False
                    break
            if create_bubble:
                bubbles.append(sphere(pos=bubble_posi, radius=bbl_r,
                        color=bubble_color,visible=True))
                bubbles[-1].popped = False
                #approximate bubble as sphere
                image_bubbles.append(sphere(pos=b_img[1]-vec(0,b_img_r,0),
                        radius=b_img_r, color=img_bubble_color,opacity=img_opac,
                        visible=True))
                image_bubbles[-1].popped = False
                
    
    for j, b in enumerate(bubbles):
        if not b.popped:
            #move the bubble up until top of radius gets to y = bubble_y_eps, easiest
            #make image and bubble invisible
            #don't have a bubble pop delay, extraneous
            #ignore change in density, radius, and speed of bubble, but not image
            b.pos = b.pos + b_v*dt
            if b.pos.y >= 0:
                b.visible = False
                b.popped = True
                image_bubbles[j].visible = False
                image_bubbles[j].popped = True
            elif b.pos.y >= -b.radius:
                image_bubbles[j].visible = False
                image_bubbles[j].popped = True
                b.visible = True
            else:   #just pop them both rather than calculate just to pop the image
                b.visible = not pov     #if pov is false, want to see, if true, don't show
                #pop bubble if it gets close to the surface
                b_impact_x = [0,0]
                b_img = [0,0]
                for i, obj_pt in enumerate([b.pos,b.pos+vec(0,b.radius,0)]):
                    [b_img[i],b_impact_x[i]] = get_image_points(obj_pt,False)                    
                
                b_img_r = mag(b_img[1]-b_img[0])
                image_bubbles[j].pos = b_img[1] - vec(0,b_img_r,0)
                image_bubbles[j].radius = b_img_r

  
    n_t += 1
    t = n_t*dt
