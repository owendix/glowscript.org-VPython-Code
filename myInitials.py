# -*- coding: utf-8 -*-
GlowScript 2.1 VPython
#Create Initials, an arc is not on glowscript but on vpython's tutorial
#I hope it is compatible

#I will make each letter of width 1, with a spacing of 0.16
#The middle initial will be centered at 0
s = 0.16
w = 1.
t = 0.37*s

#Create the O
mid=vec(-(w+s-t),0,0)
orad=0.5*w
sphere(pos=mid,radius=orad,color=vec(1,0,0))
ellipsoid(pos=mid,length=4*orad,height=orad,width=orad,
    axis=vec(0.35*w,0,1),color=vec(0,0,1))

#Create the M
dm = 0.5*w
#m1
cylinder(pos=vec(-dm+0.5*t,-dm,0),axis=vec(0,w,0),radius=t,color=vec(0.7,0.7,0))
#m2
cylinder(pos=vec(-dm,dm,0),axis=vec(dm+t,-dm-t,0),radius=t,color=vec(0,1,0))
#m3
cylinder(pos=vec(-t,-t,0),axis=vec(dm+t,dm+t,0),radius=t,color=vec(0,1,1))
#m4
cylinder(pos=vec(dm-0.5*t,dm,0),axis=vec(0,-w,0),radius=t,color=vec(0,0.5,0.1))

#Create the D
dpx = 0.5*w+s
dpy = -0.5*w
#vertical bar
cylinder(pos=vec(dpx,dpy,0),axis=vec(0,w,0),radius=t,
    color=vec(0,0,1))
#loop
#helix(pos=vec(dposx,0,0),axis=vec(0,0,0.1*t),coils=0.5,
#    radius=0.5*(w-t),thickness=2*t,up=vec(1,0,0),color=vec(1,0,1))
#dloop = curve(pos=[vec(dpx,dpy,0),vec(dpx+0.5*w,0,0),vec(dpx,-dpy,0)],radius=t,color=vec(1,0,1))
#
dl1=cylinder(pos=vec(dpx,-dpy,0),axis=vec(0.5*w,-0.5*w,0),radius=t,color=vec(1,0,1))
dl2=cylinder(pos=vec(dpx,dpy,0),axis=vec(0.5*w,0.5*w,0),radius=t,color=vec(1,0,1))