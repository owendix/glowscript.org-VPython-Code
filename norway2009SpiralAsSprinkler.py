# -*- coding: utf-8 -*-
GlowScript 2.1 VPython
"""
Models the spiral over norway in 2009 as a rotating source of expellant

Expellant stages: expellant produced, starts to evaporate, gone
Source stages: production while spinning, 
"""
#Create sprinkler 
dt = 0.025

n_arms = 3
n_red = 37
n2_red = int(n_red/3)

disp_rate=100

#how long until matter evaporates, no longer visible
t_evap = 5
evap_rate = 0.0078125   #must neatly go into 1 (picking 1/2^3 for roundoff)
#How long until black hole forms: when matter stops being ejected
t_clear = 30

#rotates clockwise
v0 = 100        #m/s (a guess at expellant speed)
Omega0 = -2.5   #rad/s, a guess at rotation rate

ang_i = 0       #initial radians
t = 0           #s
ang = ang_i
n = 0

#N particles
t_tot = t_clear + t_evap + (dt/evap_rate)
N_tot = t_tot/dt   #
the_range = (t_evap+dt/evap_rate)*v0
particle = [0 for i in range(N_tot)]

scene = display(forward=vec(0,0,-1),center=vec(0,0,0),
    range=the_range,autoscale=0)


dtheta = 0.0025
#t = 0 here
L = 400
dOmega0 = dt*Omega0
dr = v0*dt

g = graph(align="right",width=L,height=L)
yhandle=series(graph=g,label="r=a+b*theta, a=v0*t, b=-v0/Omega0")
rv = v0*t
dr_rad=dtheta*v0/Omega0
r = rv
for theta in arange(0,n_arms*2*pi,dtheta):
    r += dr_rad
    x=r*cos(theta)
    y=r*sin(theta)
    
    yhandle.plot(x,y)
        
print('https://www.youtube.com/watch?v=Ra7FMnpWMhY - start at t ~ 27s')
print('Giant sky spiral shocks Norwegians in December 2009, by AS N')
print('For comparison, this simulation produced by assuming a source of expellant')
print('which travels outward at a constant velocity v0 =',v0,'m/s. The source')
print('rotates at a constant \Omega0 =',Omega0,'rad/s, making an archimedian spiral')

while t < t_tot:
    rate(disp_rate)
    #Create matter particles (eject)
    if t < t_clear:
        s = sphere(pos=vec(0,0,0),radius=10,opacity=1.,color=color.white)
        s.v = v0*vec(cos(ang),sin(ang),0)
        if n%n_red == 13:
            s.color=color.red
        particle[n]=s
        
        #particle n, starting from 0 produced at time n*dt
        #particle (current) lifetime = t - n*dt

    i = 0
    for p in particle[0:n+1]:
        
        #Move matter particles
        #Simple forward euler numerical integration
        p.pos = p.pos + p.v*dt
        
        #Evaporate matter particles
        if t - i*dt >= t_evap and p.opacity > 0:
            p.opacity -= evap_rate    #evap_rate must be go into 1 integer times
        elif p.opacity <= 0 and p.visible == True:
            p.visible = False
        i += 1

    t += dt
    if t < t_clear:
        n += 1
    ang += dOmega0    #slight change in source angle
    
    
#Create plot: archimedian spiral:
#r = a + b*ang
#x = r*cos(ang)=(a+b*ang)*cos(ang)
#y = r*sin(ang)=(a+b*ang)*cos(ang)
#b = -v0*Omega0
#r = v0*(t - ang/Omega0)


#Annotate
#youtube videos, two
#Equations
