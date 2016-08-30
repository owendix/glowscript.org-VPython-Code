GlowScript 2.1 VPython

print("""This demonstrates the precession and nutation of a bike wheel
    (shown with a square tire for visual effect) under the influence of gravity.
    After long times, numerical error builds up and causes energy to be poorly
    conserved.
    """)


dw123 = vec(1,0,0)      #best (0,0,0) (1,0,0) (-1,0,0) for spin=50
#dw123 = vec(-5,0,0)      #initial jolt in 1,2,3 coordinates
#spin is along 3-axis, so don't choose that one

dt = 2**-8      #pick very small dt: dt = 2**-8 seems best
the_rate = 100
#t_stop = 100   #s
n_plot = 10

spin = 50        #rad/s, primary spin along 3-axis (body-fixed)

canv_width=1000

L = 0.2         #m, length of one handle/axis: half total length
g = 9.81        #m/s^2

r = 0.015       #m, radius of handle
m = 0.5         #kg, mass of handle

#Moments assume a cylinder, but a square is used for aesthetics and
#ease of vector definitions. It has the dynamics of a bike wheel
Iaxle = m*(3*r**2 + 4*L**2)/12
Iaxle3 = 0.5*m*r**2
mw = 4.5        #kg, mass of wheel/disk (consider uniform)
R = 0.3        #m, radius of wheel/disk (consider uniform)
d = 0.04        #m, thickness of wheel
Iwheel = mw*(3*R**2 + d**2)/12
Iwheel3 = 0.5*mw*R**2

#These are the body-fixed principal moments, Is around the 3-axis (spin axis),
#The other two are symmetric. The eigenvectors of this matrix are orthogonal
#(principal axis theorem), all other matrix components are zero.
I = Iwheel + Iaxle + (m+mw)*L**2
Is = Iaxle3 + Iwheel3


#aim camera along -x axis, z-axis is up
scene = display(width=canv_width,height=(9/16)*canv_width,forward=vec(-1,0,0),center=vec(0,0,0),up=vec(0,0,1),
                background=color.gray(0.7))

##Initialize alpha, beta, gamma
alpha0 = radians(-40)        #These are chosen for aesthetics, with camera = -x
beta0 = radians(35)
gamma0 = 0                  #gamma will be the primary spin axis
#Initialize
alpha = alpha0
beta = beta0
gamma = gamma0


def convert123ToXYZ(the_vec,alph,bet,gam):
    #the_vec must be a vec(#,#,#) in 123 components
    #alph,bet,gam are euler angles
    #returns a vector. I'm not sure if this function is needed
    #Checked and passed
    ca = cos(alph)
    cb = cos(bet)
    cg = cos(gam)
    sa = sin(alph)
    sb = sin(bet)
    sg = sin(gam)
    cbcg = cb*cg
    cbsg = cb*sg
    #has been checked
    Arow0 = vec(ca*cg - sa*cbsg, -(ca*sg + sa*cbcg), sa*sb)
    Arow1 = vec(sa*cg + ca*cbsg, -(sa*sg - ca*cbcg), -ca*sb)
    Arow2 = vec(sb*sg, sb*cg, cb)
    
    return vec(dot(Arow0,the_vec),dot(Arow1,the_vec),dot(Arow2,the_vec))


def getEulerAngles(ax1,ax3):
    #ax3 must be the 3-axis in xyz coordinates
    #atan2(y,x) domain: [-pi,pi], including endpoints!
    #alpha: [-pi,pi]
    #beta: [0,pi]
    #gamma: [-pi,pi]
    bet = acos(ax3.z/mag(ax3))
    if bet != 0 and bet != pi:
        xy3 = ax3 - vec(0,0,ax3.z)      #subtract projection onto z
        alph = atan2(xy3.x,-xy3.y)     #signs critical
        xp = cross(ax3,xy3)
        xp = xp/mag(xp)                 #xprime, x', unit vector
        yp = cross(ax3,xp)
        yp = yp/mag(yp)
        gam = atan2(dot(ax1,yp),dot(ax1,xp))
    else:   #know 1,2 in xy plane, use gamma = 0, alpha to x axis
        gam = 0
        alph = atan2(ax1.y,ax1.x)
    
    return alph, bet, gam

ax1 = vec(1,0,0)
ax3 = vec(0,0,1)
ax1 = convert123ToXYZ(ax1,alpha,beta,gamma)
ax3 = convert123ToXYZ(ax3,alpha,beta,gamma)

stand = cylinder(pos=vec(0,0,0),axis=vec(0,0,-3*L),radius=r,color=color.cyan)
xax = cylinder(pos=vec(0,0,-3*L),axis=vec(3*L,0,0),radius=r,color=color.magenta)
yax = cylinder(pos=vec(0,0,-3*L),axis=vec(0,3*L,0),radius=r,color=color.yellow)
#Create axle, handle end at origin
axle = cylinder(pos=vec(0,0,0),axis=2*L*ax3,radius=r,color=color.red)
#Create wheel, centerd at middle of handle
tire = box(pos=(L-0.5*d)*ax3,axis=d*ax3,length=d,width=R,height=R,up=ax1,
           color=vec(0.5,0.5,1))
#Call the whole thing the wheel
#accessing the components of the wheel is not working
#wheel=compound([axle,tire])
#print(wheel.axle.axis)         #Doesn't work, need to rotate objects individ.

w123 = dw123 + vec(0,0,spin)

t = 0   #Don't need t (autonomous ode). If I  solve analytically, I will
n = 0

IsI = (Is-I)/I
M = m + mw      #Total mass
TgI = M*g*L/I     #Ang-accel maximum
#remains constant
w3 = w123.z
#Graphs
half_canv = canv_width/2
w123_g = graph(width=half_canv,align="left")
wxyz_g = graph(width=half_canv,align="right")
L123_g = graph(width=half_canv,align="left")
Lxyz_g = graph(width=half_canv,align="right")
L2_g = graph(width=half_canv,align="left")
E_g = graph(width=half_canv,align="right")
eulerAngs_g = graph(width=canv_width,align="left")
up_g = graph(width=half_canv,align="left")
ax_g = graph(width=half_canv,align="right")
w1_s = series(graph=w123_g,label="w1 (rad/s)",color=color.red)
w2_s = series(graph=w123_g,label="w2 (rad/s)",color=color.green)
w3_s = series(graph=w123_g,label="w3 (rad/s)",color=color.blue)
wx_s = series(graph=wxyz_g,label="wx (rad/s)",color=color.red)
wy_s = series(graph=wxyz_g,label="wy (rad/s)",color=color.green)
wz_s = series(graph=wxyz_g,label="wz (rad/s)",color=color.blue)
L123x_s = series(graph=L123_g,label="L123x (kgm^2/s)",color=color.red)
L123y_s = series(graph=L123_g,label="L123y (kgm^2/s)",color=color.green)
L123z_s = series(graph=L123_g,label="L123z (kgm^2/s)",color=color.blue)
Lx_s = series(graph=Lxyz_g,label="Lx (kgm^2/s)",color=color.cyan)
Ly_s = series(graph=Lxyz_g,label="Ly (kgm^2/s)",color=color.yellow)
Lz_s = series(graph=Lxyz_g,label="Lz (kgm^2/s)",color=color.magenta)
L123_2_s = gdots(graph=L2_g,label="|L123|^2 (kg^2m^4/s^2)",color=color.red)
Lxyz_2_s = series(graph=L2_g,label="|Lxyz|^2 (kg^2m^4/s^2)",color=color.green)
Ek123_s = gdots(graph=E_g,label="Ek_123 (J)",color=color.blue)
Ekxyz_s = series(graph=E_g,label="Ek_xyz (J)",color=color.magenta)
alph_s = series(graph=eulerAngs_g,label="alpha (deg)",color=color.red)
bet_s = series(graph=eulerAngs_g,label="beta (deg)",color=color.green)
gam_s = series(graph=eulerAngs_g,label="gamma (deg)",color=color.blue)
upx_s = series(graph=up_g,label="up_x (m)",color=color.red)
upy_s = series(graph=up_g,label="up_y (m)",color=color.green)
upz_s = series(graph=up_g,label="up_z (m)",color=color.blue)
axx_s = series(graph=ax_g,label="ax_x (m)",color=color.red)
axy_s = series(graph=ax_g,label="ax_y (m)",color=color.green)
axz_s = series(graph=ax_g,label="ax_z (m)",color=color.blue)


while True:
    #Numerically integrate: midpoint. Simultaneously updates angles, omegas
    rate(the_rate)
    
    #Update angles
    ang = mag(w123)*2*dt
    wxyz = convert123ToXYZ(w123,alpha,beta,gamma)
    #wxyz = dalpha/dt + dbeta/dt + dgamma/dt
    #graph
    if n%n_plot == 0:
        L123 = vec(I*w123.x,I*w123.y,Is*w123.z)
        Lxyz = convert123ToXYZ(L123,alpha,beta,gamma)
        Ek123 = 0.5*dot(w123,L123)
        Ekxyz = 0.5*dot(wxyz,Lxyz)
        w1_s.plot(t,w123.x)
        w2_s.plot(t,w123.y)
        w3_s.plot(t,w123.z)
        wx_s.plot(t,wxyz.x)
        wy_s.plot(t,wxyz.y)
        wz_s.plot(t,wxyz.z)
        L123x_s.plot(t,L123.x)
        L123y_s.plot(t,L123.y)
        L123z_s.plot(t,L123.z)
        Lx_s.plot(t,Lxyz.x)
        Ly_s.plot(t,Lxyz.y)
        Lz_s.plot(t,Lxyz.z)
        L123_2_s.plot(t,mag2(L123))
        Lxyz_2_s.plot(t,mag2(Lxyz))
        Ek123_s.plot(t,Ek123)
        Ekxyz_s.plot(t,Ekxyz)
        alph_s.plot(t,degrees(alpha))
        bet_s.plot(t,degrees(beta))
        gam_s.plot(t,degrees(gamma))
        upx_s.plot(t,tire.up.x)
        upy_s.plot(t,tire.up.y)
        upz_s.plot(t,tire.up.z)
        axx_s.plot(t,axle.axis.x)
        axy_s.plot(t,axle.axis.y)
        axz_s.plot(t,axle.axis.z)

    axle.rotate(angle=ang,axis=hat(wxyz),origin=vec(0,0,0))
    tire.rotate(angle=ang,axis=hat(wxyz),origin=vec(0,0,0))

    alpha, beta, gamma = getEulerAngles(tire.up,axle.axis)


#update omega
sb = sin(beta)
    cg = cos(gamma)
    sg = sin(gamma)
    TgIsb = TgI*sb
    
    w1 = w123.x
    w2 = w123.y
    if n == 0:
        w1old = w1
        w2old = w2
        dt = 0.5*dt
"""
    #Torque free
    tmp1 = w1old - w2*w3*IsI*2*dt
    tmp2 = w2old + w3*w1*IsI*2*dt
    #w3 is unchanged
    """
        #Torqued by gravity
        tmp1 = w1old + (TgIsb*cg - w2*w3*IsI)*2*dt
            tmp2 = w2old + (-TgIsb*sg + w3*w1*IsI)*2*dt
#w3 is unchanged

w1old = w1
    w2old = w2
    
    #update book axis and up: use rotate? or axis? Try rotate
    w1 = tmp1
    w2 = tmp2
    if n == 0:
        dt = 2*dt
#now update angle (mixes future omega with current angle)
#reupdate w123 for convenience
w123 = vec(w1,w2,w3)
    
    n += 1
    t += dt
