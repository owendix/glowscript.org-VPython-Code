GlowScript 2.1 VPython

print("""
    Simulate free rotation of a uniform rectangular solid w dimensions a > b > c
    This demonstrates that spinning around the intermediate principal axis (b)
    is unstable.
    """)
"""
    This solves with both the analytical and numerical methods but the
    analytical solution is a linearized approximation in the vicinity of the small
    wobbling initial conditions, so its long term behavior is not trustworthy.
    In both cases,
    the extreme axes are stable and the intermediate axis is unstable. Both,
    however, show unphysical instability when using this axis, though in different
    ways. I have tested the components of this algorithm, including plotting
    the euler angles and various other vector components.
    Although there is always a possibility of an error, there would need to be one
    in both the numerical and analytical simulation. I believe the odd
    unstable behavior for the numeric simulation arrives from the unstable
    numerical integration technique. The odd unstable behavior with the
    analytical method is because it is linearized in the vicinity of the initial
    conditions so its long term behavior is not to be trusted. The exponential that
    appears is not tapered off so the behavior we see is exactly what we should
    expect from such an approximation.
    
    It would be quite onerous to implement a stable numerical algorithm, such
    as the implicit midpoint, because this involves optimization such as with
    Newton's method, and this is much easier with linear algebra technique like
    those found in NumPy or SciPy...somethings glowscript does not let you import.
    It is possible that this could be fixed by forcing angular momentum or energy
    in the inertial reference frame to be constant, but the dynamics are good
    enough for me. They demonstrate the instability and correct use of
    Euler angles and the Euler equation. The numerical solution
    uses a midpoint type numerical integration algorithm with choices in
    how to update the different parameters. The best theoretically uses
    simultaneous updates. (which_method = 2).
    """

#Initial conditions
spin_axis = ['a','c','b']       #'a', 'b', or 'c': primary spin axis, all perturbed

solveNumerically = False
print('Solving analytically? ',not solveNumerically)

canv_width=1000

dt = 2**-5      #pick very small dt: dt = 2**-4 --> no tumble any axis, 2**-3 all unstable
the_rate = 1000
t_stop = 100   #s
n_plot = 10

spin = 1        #rad/s, primary spin along one spin_axis
dw = 0.05       #rad/s, perturbation

#Dimensions of the uniform rectangular solid
#must have a > b > c
a = 0.236       #m
b = 0.181       #m
c = 0.029       #m

m = 0.804       #kg

scene = display(width=canv_width,height=(9/16)*canv_width,forward=vec(-1,0,0),center=vec(0,0,0),up=vec(0,0,1),
                background=color.gray(0.7))

w123_g = graph(width=500,align="left")
wxyz_g = graph(width=500,align="right")
eulerAngs_g = graph(width=canv_width,align="left")
up_g = graph(width=500,align="left")
ax_g = graph(width=500,align="right")
w1_s = series(graph=w123_g,label="w1",color=color.red)
w2_s = series(graph=w123_g,label="w2",color=color.green)
w3_s = series(graph=w123_g,label="w3",color=color.blue)
wx_s = series(graph=wxyz_g,label="wx",color=color.red)
wy_s = series(graph=wxyz_g,label="wy",color=color.green)
wz_s = series(graph=wxyz_g,label="wz",color=color.blue)
alph_s = series(graph=eulerAngs_g,label="alpha",color=color.red)
bet_s = series(graph=eulerAngs_g,label="beta",color=color.green)
gam_s = series(graph=eulerAngs_g,label="gamma",color=color.blue)
upx_s = series(graph=up_g,label="up_x",color=color.red)
upy_s = series(graph=up_g,label="up_y",color=color.green)
upz_s = series(graph=up_g,label="up_z",color=color.blue)
axx_s = series(graph=ax_g,label="ax_x",color=color.red)
axy_s = series(graph=ax_g,label="ax_y",color=color.green)
axz_s = series(graph=ax_g,label="ax_z",color=color.blue)

myseries=[w1_s,w2_s,w3_s,wx_s,wy_s,wz_s,alph_s,bet_s,gam_s,upx_s,upy_s,upz_s,
          axx_s,axy_s,axz_s]

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
scount = 0

for s_ax in spin_axis:
    scount += 1
    print('Primary spin axis aligned with: '+spin_axis+'-dimension (a>b>c)')
    print('Other axes are perturbed')
    #aim camera along -x axis, z-axis is up
    for srs in myseries:
        srs.delete()
    
    ##Initialize alpha, beta, gamma
    alpha0 = radians(-40)        #These are chosen for aesthetics, with camera = -x
    beta0 = radians(35)
    gamma0 = 0                  #gamma will be the primary spin axis
    #Initialize
    alpha = alpha0
    beta = beta0
    gamma = gamma0
    
    ax1 = vec(1,0,0)
    ax2 = vec(0,1,0)
    ax3 = vec(0,0,1)
    ax1 = convert123ToXYZ(ax1,alpha,beta,gamma)
    ax2 = convert123ToXYZ(ax2,alpha,beta,gamma)
    ax3 = convert123ToXYZ(ax3,alpha,beta,gamma)
    
    #All other matrix components are zero. Because these are fixed to the body,
    #these principal moments will not change with time
    #The fact that they are principal moments means the principal axes are perpendicular
    #Principal axis theorem: the rotational inertia matrix will be orthogonal
    
    #These assignments are cylical (a,b,c) and (1,2,3)
    #I checked, visually. They work...so far so good I think
    #Initialize orientation alpha,beta, gamma
    if scount > 1:
        book.visible=False
    if 'a' in s_ax:
        print('a -->3, b -->1, c-->2')
        I3 = m*(b**2 + c**2)/12        #1-vec points along largest, a, axis
        I1 = m*(a**2 + c**2)/12        #2-vec points along intermediate, b, axis
        I2 = m*(a**2 + b**2)/12        #3-vec points along shortest, c, axis
        book = box(pos=vec(0,0,0),axis=ax3,length=a,height=b,width=c,up=ax1)
        the_ln, the_ht, the_wd = a, b, c
    elif 'b' in s_ax:
        print('b-->3, c-->1, a-->2')
        I2 = m*(b**2 + c**2)/12        #1-vec points along largest, a, axis
        I3 = m*(a**2 + c**2)/12        #2-vec points along intermediate, b, axis
        I1 = m*(a**2 + b**2)/12        #3-vec points along shortest, c, axis
        book = box(pos=vec(0,0,0),axis=ax3,length=b,height=c,width=a,up=ax1)
        the_ln, the_ht, the_wd = b, c, a
    else:   #'c' in spin axis
        print('c-->3, a-->1, b-->2')
        #Fixed-body, center of mass origin, principal axes-basis moments for a
        #uniform rectangular solid with dimensions a > b > c
        I1 = m*(b**2 + c**2)/12        #1-vec points along largest, a, axis
        I2 = m*(a**2 + c**2)/12        #2-vec points along intermediate, b, axis
        I3 = m*(a**2 + b**2)/12        #3-vec points along shortest, c, axis
        book = box(pos=vec(0,0,0),axis=ax3,length=c,height=a,width=b,up=ax1)
        the_ln, the_ht, the_wd = c, a, b
    #w = w + vec(dw,dw,dw)   #perturb the spin a little bit, need in x,y,z basis
    book.color=vec(0.5,0.5,1)   #light blue
    w123 = vec(dw,dw,spin+dw)   #initialize omega in 123 coordinates, spin=3-axis

    #establish vectors, visually: w = omega = d(ang)/dt
    #Need to learn how VPython rotates things
    #vectors have x,y,z components but need to rotate this to the
    #basis attached to the rectangular solid
    #may need to use Euler angles
    #or just maintain the 1,2,3 basis separately???

    #book.rotate(angle=a_in_rads, axis=vec(x,y,z))
    t = 0   #Don't need t (autonomous ode). If I  solve analytically, I will
    n = 0

    I321 = (I3-I2)/I1
    I132 = (I1-I3)/I2
    I213 = (I2-I1)/I3
    Ifrac1 = I132*I321

    if solveNumerically:
        while t < t_stop:
            #Numerically integrate: midpoint. Simultaneously updates angles, omegas
            rate(the_rate)
            
            ang = mag(w123)*2*dt
            wxyz = convert123ToXYZ(w123,alpha,beta,gamma)
            
            if n%n_plot == 0:
                w1_s.plot(t,w123.x)
                w2_s.plot(t,w123.y)
                w3_s.plot(t,w123.z)
                wx_s.plot(t,wxyz.x)
                wy_s.plot(t,wxyz.y)
                wz_s.plot(t,wxyz.z)
                alph_s.plot(t,degrees(alpha))
                bet_s.plot(t,degrees(beta))
                gam_s.plot(t,degrees(gamma))
                upx_s.plot(t,book.up.x)
                upy_s.plot(t,book.up.y)
                upz_s.plot(t,book.up.z)
                axx_s.plot(t,book.axis.x)
                axy_s.plot(t,book.axis.y)
                axz_s.plot(t,book.axis.z)
        
            #wxyz = dalpha/dt + dbeta/dt + dgamma/dt
            book.rotate(angle=ang,axis=hat(wxyz))
            #update the axis and up vectors?
            
            #update angles
            alpha, beta, gamma = getEulerAngles(book.up,book.axis)
            
            w1 = w123.x
            w2 = w123.y
            w3 = w123.z
            if n == 0:
                w1old = w1
                w2old = w2
                w3old = w3
                dt = 0.5*dt

    tmp1 = w1old-w2*w3*I321*2*dt
        tmp2 = w2old-w3*w1*I132*2*dt
            tmp3 = w3old-w1*w2*I213*2*dt
            
            w1old = w1
            w2old = w2
            w3old = w3
            
            #update book axis and up: use rotate? or axis? Try rotate
            w1 = tmp1
            w2 = tmp2
            w3 = tmp3
            if n == 0:
                dt = 2*dt
        #now update angle (mixes future omega with current angle)
        #reupdate w123 for convenience
        w123 = vec(w1,w2,w3)
            
            n += 1
            t += dt

else:
    #analytical solution, obtained by assuming a small perturbation
    #along the 1 and 2 direction (can be different magnitudes)
    #then linearizing, turning into a matrix-differential equation and solving
    #w123 contains the initial angular velocities in the 123-basis
    w1 = w123.x
        w2 = w123.y
        w3 = w123.z #constant: this lets you linearize
        
        #(vec) w(t) = C_p*u_p*exp(lam_p*t) + C_m*u_m*exp(lam_m*t)
        #u_p, u_m are eigenvectors, lam_p, lam_m are eigenvalues
        #C_p, C_m are constants determined by initial conditions
        #see above, naming convention I132 = (I1 - I3)/I2 in that order
        
        
        t = 0
        if Ifrac1 > 0:   #happens when I3 is intermediate
            lam_p = abs(w3)*sqrt(Ifrac1)
            lam_m = -lam_p
            u_p = -w3*I321/lam_p
            u_m = -u_p
            #sqIfrac2 = sqrt(Ifrac2)
            tmp = w1/u_m
            C_p = -0.5*(tmp - w2)    #only want to use the initial values
            C_m = 0.5*(tmp + w2)
            Cpup = C_p*u_p
            Cmum = C_m*u_m
            
            while t < t_stop:
                rate(the_rate)
                e_p = exp(lam_p*t)
                e_m = exp(lam_m*t)
                w1 = Cpup*e_p + Cmum*e_m   #I suppose when this blows up, linear approx is bad
                w2 = C_p*e_p + C_m*e_m
                
                alpha, beta, gamma = getEulerAngles(book.up,book.axis)
                w123=vec(w1,w2,w3)
                ang = mag(w123)*dt
                wxyz = convert123ToXYZ(w123,alpha,beta,gamma)
                
                if n%n_plot == 0:
                    w1_s.plot(t,w123.x)
                    w2_s.plot(t,w123.y)
                    w3_s.plot(t,w123.z)
                    wx_s.plot(t,wxyz.x)
                    wy_s.plot(t,wxyz.y)
                    wz_s.plot(t,wxyz.z)
                    alph_s.plot(t,degrees(alpha))
                    bet_s.plot(t,degrees(beta))
                    gam_s.plot(t,degrees(gamma))
                    upx_s.plot(t,book.up.x)
                    upy_s.plot(t,book.up.y)
                    upz_s.plot(t,book.up.z)
                    axx_s.plot(t,book.axis.x)
                    axy_s.plot(t,book.axis.y)
                    axz_s.plot(t,book.axis.z)
                
                book.rotate(angle=ang,axis=hat(wxyz))
                
                t += dt

    else:           #happens when I3 is an extreme
        #both fractions imaginary...exp(i*lam_p*t) = vec(cos,sin,0), take Real part
        lam_p = abs(w3)*sqrt(abs(Ifrac1))
            u_m = w3*I321/lam_p
            #sqIfrac2 = sqrt(Ifrac2)
            w1cos = 0.5*w1
            w1sin = -0.5*w2*u_m
            w2cos = 0.5*w2
            w2sin = 0.5*w1/u_m
            while t < t_stop:
                rate(the_rate)
                
                #w1 = C_p*sqIfrac2*sin(lam_p*t) - C_m*sqIfrac2*sin(lam_m*t)
                #w2 = C_p*cos(lam_p*t) + C_m*cos(lam_m*t)
                w1 = 2*(w1cos*cos(lam_p*t) + w1sin*sin(lam_p*t))
                w2 = 2*(w2cos*cos(lam_p*t) + w2sin*sin(lam_p*t))
                
                alpha, beta, gamma = getEulerAngles(book.up,book.axis)
                w123=vec(w1,w2,w3)
                ang = mag(w123)*dt
                wxyz = convert123ToXYZ(w123,alpha,beta,gamma)
                
                if n%n_plot == 0:
                    w1_s.plot(t,w123.x)
                    w2_s.plot(t,w123.y)
                    w3_s.plot(t,w123.z)
                    wx_s.plot(t,wxyz.x)
                    wy_s.plot(t,wxyz.y)
                    wz_s.plot(t,wxyz.z)
                    alph_s.plot(t,degrees(alpha))
                    bet_s.plot(t,degrees(beta))
                    gam_s.plot(t,degrees(gamma))
                    upx_s.plot(t,book.up.x)
                    upy_s.plot(t,book.up.y)
                    upz_s.plot(t,book.up.z)
                    axx_s.plot(t,book.axis.x)
                    axy_s.plot(t,book.axis.y)
                    axz_s.plot(t,book.axis.z)
                
                book.rotate(angle=ang,axis=hat(wxyz))
                
                t += dt