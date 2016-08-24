# -*- coding: utf-8 -*-
GlowScript 2.1 VPython

dt = 200
the_rate = 2000    #2000 good for viewing

#Ignores asteroid presence: sets asteroid.mass=1e-7, 
#asteroid.color=color.black, asteroid.opacity=0
ignoreAsteroid = True   #True, False


#moon properties
smaj = 3.844e8      #meters, semimajor axis, starts at positive y
vbase = 0.7e3
#We want apogee: semimajor + distance to focus
apogee = 4.055e8      #meters, apogee, farthest from earth
vmin = 964        #m/s, minimum speed


rinit = apogee      #apogee or smaj
if rinit == apogee:
    vinit = vmin
    erasefact = 24
else:
    rinit = smaj
    vinit = vbase
    erasefact = 14

meth_idx = 0    #0, 1, 2, 3
#Checked with asteroid off (ignoreAsteroid=True)
#0, 2: both pretty stable but display error (black or white regions) appear
#1: simultaneous euler is not good for conserving energy, mixed is better
#3: stable but display error occurs, maybe as fast or faster than 0, 2

methods = ['mixedeuler','simulteuler','midpt','rk']
method = methods[meth_idx]


usePlot = True  #True, False
nVecDraw = 500
nVecErase = erasefact*nVecDraw     #14, or 24
nPlot = 100


#things to graph vs time: 
#swept area, total energy, kinetic energy, gravitational energy
gArea = graph(ymin=0)
gEnergy = graph()
swept_area = series(graph=gArea,color=color.magenta,label='Swept Area(m^2)')
Etot = series(graph=gEnergy,color=color.red,label='E_total(J)')
Ekin = series(graph=gEnergy,color=color.blue,label='Ek(J)')
Eg = series(graph=gEnergy,color=color.green,label='Eg(J)')

def drawmotionvecs(obj,origin=earth.pos):
    #These objects have methods: pos, p
    #Could later add accel. So far this seems easier than
    #the build in make_arrow feature, which has little
    #documentation
    C = 0.005    #scale for shaftwidth
    rmag = mag(obj.pos)
    sw = C*rmag
    #rpscale = 1e-15
    rvec = arrow(pos=origin,axis=obj.pos,shaftwidth=sw)
    #pvec = arrow(pos=obj.pos,axis=rpscale*hat(obj.p),shaftwidth=sw)
    #Momentum has too widely varying scale compared to r to map both and 
    #have significance in the momentum magnitude
    #return [rvec, pvec]
    
    return rvec
    
def addLsVecs(a,b):
    f=[]
    for i,q in enumerate(a):
        f.append(q+b[i])
    return f
    
def subLsVecs(a,b):
    f=[]
    for i,q in enumerate(a):
        f.append(q-b[i])
    return f
    
def multLsVecs(h,a):   
    return [h*q for q in a]
"""
Class operator overrides don't work with glowscript,
Copied working python class definition here:
class ListOfVecs:
    #Creates a list of vecs, for use with glowscript
    #Construct: varname = ListOfVecs(vec1,vec2,...)
    #Return list: varname.lov
    #Can add, subtract with other lists, using +, -
    #Can left and right multiply with a scalar, using *
    
    lov = [] #Access that list with varname.lv
    
    def __init__(self,*args):
        self.lov = list(args)

    def __add__(self,b):
        #Assumes self, b equal length
        s = [q+b.lov[i] for i,q in enumerate(self.lov)]
        return ListOfVecs(*s)
        
    def __sub__(self,b):
        #Assumes self, b equal length
        s = [q-b.lov[i] for i,q in enumerate(self.lov)]     
        return ListOfVecs(*s)
        
    def __mul__(self,h): 
        #Scalar h
        s=[h*q for q in self.lov]   
        return ListOfVecs(*s)
    
    def __rmul__(self,h):
        #Scalar h
        s=[h*q for q in self.lov]    
        return ListOfVecs(*s)
        
Along with (in 'mixedeuler' section):
        Y0 = ListOfVecs(asteroid.p, moon.p, asteroid.pos, moon.pos)
        F = ListOfVecs(FtotalA,FtotalM,
                Y0.lov[0]/asteroid.mass,Y0.lov[1]/asteroid.mass)
        Y1 = Y0 + F*dt
        [asteroid.p,moon.p,asteroid.pos,moon.pos]=Y1.lov
        
Result --> Error: Y0.+ is not a function
"""

def Forces(apos=asteroid.pos, mpos=moon.pos, epos=earth.pos):
    #EM calc
    rEM=epos-mpos
    rEMmag=mag(rEM)
    FEonM=(G*moon.mass*earth.mass/(rEMmag**3))*rEM

    #EA calc
    rEA=epos-apos
    rEAmag=mag(rEA)
    FEonA=(G*asteroid.mass*earth.mass/(rEAmag**3))*rEA

    #MA calc
    rMA=mpos-apos
    rMAmag=mag(rMA)
    FMonA=(G*asteroid.mass*moon.mass/(rMAmag**3))*rMA

    #Ftotal moon
    FtotalM=FEonM-FMonA

    #Fotal asteroid
    FtotalA=FEonA+FMonA

    return [FtotalA, FtotalM]

def Egtotal(apos=asteroid.pos, mpos=moon.pos, epos=earth.pos):
    ##EM calc
    rEM=epos-mpos
    egEM=moon.mass*earth.mass/mag(rEM)

    #EA calc
    rEA=epos-apos
    egEA=asteroid.mass*earth.mass/mag(rEA)

    #MA calc
    rMA=mpos-apos
    egMA=asteroid.mass*moon.mass/mag(rMA)

    return -G*(egEM+egEA+egMA)

def sweptArea(dt,rf=moon.pos, pf=moon.p, m=moon.mass):
    """
    Approximate the incremental (dt) area previously swept in a time step
    For area yet to be swept, put before, and change first - to positive
    No need to store previous location: uses momentum, assumes dr=(p/m)*dt 
    Approximates arc as a triangle, finds area with numerically stable
    Heron's formula, and assumes |dr| is the smallest side
    """
    dr = pf*(-dt/m)
    ri = rf + dr
    
    a = mag(ri)
    b = mag(rf)
    if a < b:
        tmp = a
        a = b
        b = tmp
    c = mag(dr) #a>=b>c
    ab = a - b
    
    return 0.25*sqrt((a + (b + c))*(c - ab)*(c + ab)*(a + (b - c)))

##Below is modified from Dwain Desbien's orbit.py from PHY121

#Earth information... Assume fixed
earth = sphere(pos=vec(0,0,0), color=vec(0.7,0.7,1), radius=6.3e6)
earth.mass=6e24
#Moon Information
moon=sphere(pos=vec(0,rinit,0), color=color.white,radius=6e6, make_trail=True,
    trail_radius=1e6)
    #trail_type="points",interval=100)
moon.mass=14.6e22
moon.p=vec(vinit*moon.mass,0,0)
#satelite information
asteroid=sphere(pos=vec(0,-2e8,0),color=vec(1,0.7,0.2),radius=2e6, 
    make_trail=True,trail_radius=1e6)#,trail_type="points",interval=100)
if not ignoreAsteroid:
    asteroid.mass=8e21
else:
    asteroid.mass=1e-7
    asteroid.color=color.black
    asteroid.make_trail=False
    asteroid.visible=False
    asteroid.trail_opacity=0
    asteroid.opacity=0
    
asteroid.p=vec(-1e3*asteroid.mass,0,0)

G=6.7e-11

scene.autoscale=0


t = 0
n = 0
motionvecs = []
while True:
    rate(the_rate)
    #don't need to update time because I stop whenever I'm done
    [FtotalA, FtotalM] = Forces(asteroid.pos, moon.pos, earth.pos)

    n += 1
    #Draw position vectors periodically, shows roughly equal area, within
    #Simulation error (looks good enough)
    if n%nVecDraw == 0:
        av = drawmotionvecs(asteroid)
        av.color=asteroid.color
        mv=drawmotionvecs(moon)
        mv.color=moon.color
        motionvecs.extend([av,mv])
    if n%(nVecErase) == 0:
        for m in motionvecs:
            m.visible=False
        motionvecs = []
    
    #Update moon
    if 'mixedeuler' in method:
        moon.p=moon.p+FtotalM*dt
        moon.pos=moon.pos+(moon.p/moon.mass)*dt
        #Update asteroid
        asteroid.p = asteroid.p + FtotalA*dt
        asteroid.pos = asteroid.pos + asteroid.p*(dt/asteroid.mass)
    elif 'euler' in method:
        #Create list of vectors
        #y = [pa, pm, ra, rm]
        #Y0 = [asteroid.p, moon.p, asteroid.pos, moon.pos]
        #dY/dt = F(t,Y(t)) = [dpa/dt, dpm/dt, dra/dt, drm/dt] 
        #Effect of simultaneous updates: uncomment two lines below
        #Unstable energy, it is correctly done though
        
        Y0 = [asteroid.p, moon.p, asteroid.pos, moon.pos]
        Y1 = [Y0[0]+FtotalA*dt, Y0[1]+FtotalM*dt, 
            Y0[2]+Y0[0]*(dt/asteroid.mass), Y0[3]+Y0[1]*(dt/moon.mass)]
        
        [asteroid.p, moon.p, asteroid.pos, moon.pos] = Y1
        
    elif 'midpt' in method:
        #Effect of simultaneous updates and explicit midpoint
        #To calculate even explicit midpoint, I need to recalculate all of the forces
        Y0 = [asteroid.p, moon.p, asteroid.pos, moon.pos]
        h2 = 0.5*dt
        F0 = [FtotalA, FtotalM, Y0[0]/asteroid.mass, Y0[1]/moon.mass]
        Y1arg = addLsVecs(Y0,multLsVecs(h2,F0))
        [FtotalA, FtotalM] = Forces(Y1arg[2],Y1arg[3],earth.pos)
        F1 = [FtotalA, FtotalM, Y1arg[0]/asteroid.mass, Y1arg[1]/moon.mass]

        [asteroid.p, moon.p, asteroid.pos, moon.pos] = addLsVecs(Y0,
                                                            multLsVecs(dt,F1))

    else:
        #Effect of simultaneous updates and explicit 4th order runge-kutta
        Y0 = [asteroid.p, moon.p, asteroid.pos, moon.pos]
        F0 = [FtotalA, FtotalM, Y0[0]/asteroid.mass, Y0[1]/moon.mass]
        k1=multLsVecs(dt,F0)
        Yk = addLsVecs(Y0,multLsVecs(0.5,k1))
        [FtotalA, FtotalM] = Forces(Yk[2],Yk[3],earth.pos)
        k2 = multLsVecs(dt,[FtotalA,FtotalM,Yk[0]/asteroid.mass,Yk[1]/moon.mass])
        Yk = addLsVecs(Y0,multLsVecs(0.5,k2))
        [FtotalA, FtotalM] = Forces(Yk[2],Yk[3],earth.pos)
        k3 = multLsVecs(dt,[FtotalA,FtotalM,Yk[0]/asteroid.mass,Yk[1]/moon.mass])
        Yk = addLsVecs(Y0,k3)
        [FtotalA, FtotalM] = Forces(Yk[2],Yk[3],earth.pos)
        k4 = multLsVecs(dt,[FtotalA,FtotalM,Yk[0]/asteroid.mass,Yk[1]/moon.mass])
        Y1 = addLsVecs(Y0,multLsVecs(1/6,addLsVecs(k1,
            addLsVecs(multLsVecs(2,addLsVecs(k2,k3)),k4))))
        
        [asteroid.p, moon.p, asteroid.pos, moon.pos] = Y1
    
    if usePlot and n%nPlot == 0:
        area=sweptArea(dt)
        Ekval = 0.5*(mag2(asteroid.p)/asteroid.mass+mag2(moon.p)/moon.mass)
        Egval = Egtotal()
        
        #swept_area, Etot, Ekin, Eg
        Etot.plot(t,Ekval+Egval)
        Ekin.plot(t,Ekval)
        Eg.plot(t,Egval)
        swept_area.plot(t,area)
    
    
    t += dt
    
    #Implicit Midpoint (best, symplectic, A and L stable) 
    #is not practical with glowscript --> no numpy for 
    #matrix manipulation to solve a Jacobian to optimize implicit methods
    #Also no Scipy, which has a built in implicit integrator