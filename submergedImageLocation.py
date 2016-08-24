# -*- coding: utf-8 -*-
GlowScript 2.1 VPython
#print(asin(0.5)*180/pi)
#This finds the location of the center of an underwater object at 
#different angles. It's always closer to the light entry point, in 
#both x and y, than the object.
#Input does not work
#N = int(input("How many different angles to display? "))

N = 5
#set to 0 to leave all trails
#how often to erase the trails, best: either 0 or 1
Nerase=0
#higher number, simulates faster
therate=20000

print("This takes two nearby light rays, that could conceivably both extend from")
print("a fish and enter a person's eye, applies the law of refraction at the")
print("water/air interface, then traces virtual rays back along the line of sight")
print("to where they intersection, where a virtual image should be.")
print("\nNotice the location of the image for the object location viewed at that angle.")
print("\nKey:")
print("Blue = water, black=air")
print("Blue ellipsoid = fish's location")
print("Orange ellipsoid = fish's virtual image when viewed from angle above water")
print("Yellow/red lines are virtual rays, traced backward along straight line")
print("\nThe image is closer to the light ray entry point, in both x & y, than the object")

leftx = -2
bigy = 2

scene = display(forward=vec(0,0,-1),center=vec(0.5*leftx,0,0),
    range=-0.8*leftx) #Tilts camera

#scene.fullscreen = True
#scene.autoscale=0

n_air = 1.
n_water = 4/3
n_ratio = n_water/n_air

#Fish depth
d = -1 #meters
#Horizontal distance of second ray at water level
s = 0.001 #meters

#Draw water level: pos = middle of box

waterlevel = box(pos=vec(0.5*leftx,0.6*d,0),
    length=1.2*abs(leftx),height=-1.2*d,width=3,color=color.blue,opacity=0.2)
    
#Draw fish object, light blue
l2 = 0.15
fish = ellipsoid(pos=vec(0,d,0),length=2*l2,height=0.1,width=0.1,
    color=color.blue)
#tail = triangle(v0=vertex(pos=vec(0,d,0),color=color.blue),
#        v1=vertex(pos=vec(0.3,d+0.05,0),color=color.blue),
#        v2=vertex(pos=vec(0.3,d-0.05,0),color=color.blue))
#fish = compound([body,tail])
#objtext = text(text='O',pos=vec(0,y,0),align='center',length=0.05,
#    height=0.05,color=color.green)
ray1=sphere(pos=fish.pos,radius=0.1*fish.width,color=color.white,
    make_trail=True,trail_radius=0.1*fish.width,opacity=0.2)#retain=10)
ray2=sphere(pos=fish.pos,radius=0.1*fish.width,color=color.green,
    make_trail=True,trail_radius=0.1*fish.width,opacity=0.2)#,retain=10)

#theta range in radians (0.01745rad = 1deg)
ang_i = 0.1
#critical angle in radians
ang_crit = asin(1/n_ratio)
dang = (ang_crit-2*ang_i)/N

#scale dang for color scheme
dcolor = dang/(ang_crit-ang_i)

ang = ang_i
n = 0

while ang < (ang_crit-dang):
    n+=1
    
    rate(therate)

    xd = -d*tan(ang)    #positive, d neg, fish at (0,d,0)
    
    #calculate angles
    ang_s = atan((s+xd)/-d)
    ang_img = asin(n_ratio*sin(ang))
    ang_s_img = asin(n_ratio*sin(ang_s))
    #calculate image location
    m = -1/tan(ang_img)
    m_s = -1/tan(ang_s_img)  
    unshiftx = m_s*s/(m_s - m)
    imgx = -xd-unshiftx
    imgy = -m*unshiftx


    #Show light rays, travel at water speeds
    ray1.v = vec(-xd,-d,0)
    ray1.v = hat(ray1.v)/n_water
    #second ray hits water just a bit left of first
    ray2.v = vec(-xd-s,-d,0)
    ray2.v = hat(ray2.v)/n_water
    
    #step size for light
    dt = 0.0001
    
    #Light is under the water the water
    while ray2.pos.y < 0:
        rate(therate)
        if ray1.pos.y < 0:
            ray1.pos = ray1.pos+ray1.v*dt
        ray2.pos = ray2.pos+ray2.v*dt
    #Breaks when it is greater than 0
    
    #Extrapolate to the water (need to be exactly there)
    ray1.pos = ray1.pos - ray1.v*(ray1.pos.y/ray1.v.y)
    ray2.pos = ray2.pos - ray2.v*(ray2.pos.y/ray2.v.y)
    
    #Light has left the water
    ray1.v = vec(-sin(ang_img),cos(ang_img),0)/n_air
    ray2.v = vec(-sin(ang_s_img),cos(ang_s_img),0)/n_air
    
    
    #Create virtual image light rays, once
    if (ang < ang_i+0.5*dang):
        virt1=sphere(pos=ray1.pos,radius=0.1*fish.width,color=color.yellow,
            make_trail=True,trail_radius=0.1*fish.width,opacity=0.2)
        virt2=sphere(pos=ray2.pos,radius=0.1*fish.width,color=color.red,
            make_trail=True,trail_radius=0.1*fish.width,opacity=0.2)

    virt1.make_trail=False
    virt2.make_trail=False
    virt1.pos=ray1.pos
    virt2.pos=ray2.pos
    virt1.make_trail=True
    virt2.make_trail=True
    virt1.v=-ray1.v
    virt2.v=-ray2.v

    #Light is above the water
    while ray1.pos.x > leftx and ray2.pos.y < bigy:
        rate(therate)
        ray1.pos = ray1.pos+ray1.v*dt
        ray2.pos = ray2.pos+ray2.v*dt
    #Virtual image light rays extend below
    while virt1.pos.y > imgy or virt2.pos.y > imgy:
        rate(therate)
        if virt1.pos.y > imgy:
            virt1.pos = virt1.pos+virt1.v*dt
        if virt2.pos.y > imgy:
            virt2.pos = virt2.pos+virt2.v*dt
    
    if ang < ang_i+0.5*dang:
        #Create image once and move, not multiple times
        image = ellipsoid(pos=vec(imgx,imgy,0),length=2*l2,height=0.1,width=0.1,
            color=color.orange,opacity=0.9)
        #imagetail = triangle(v0=vertex(pos=vec(imgx,imgy,0),color=color.orange,
        #        opacity=0.5),
        #        v1=vertex(pos=vec(imgx+0.3,imgy+0.05,0),color=color.orange,
        #        opacity=0.5),
        #        v2=vertex(pos=vec(imgx+0.3,imgy-0.05,0),color=color.orange,
        #        opacity=0.5))
        #image = compound([imagebody,imagetail])
    #Move image of fish
    image.pos=vec(imgx,imgy,0)
           
    ray1.make_trail=False
    ray2.make_trail=False
    ray1.pos=fish.pos
    ray2.pos=fish.pos
    
    #Erase every so often: don't this erases everything
    if n%Nerase==0:
        ray1.clear_trail()
        ray2.clear_trail()
        virt1.clear_trail()
        virt2.clear_trail()    
    ray1.make_trail=True
    ray2.make_trail=True
    #ray1.color=vec(1,ray1.color.y-dcolor,ray1.color.z-dcolor)
    #ray2.color=vec(ray2.color.x-dcolor,ray2.color.y-dcolor,1)
    ang += dang