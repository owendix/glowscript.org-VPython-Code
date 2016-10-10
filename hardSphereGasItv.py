GlowScript 2.1 VPython
#Ideal gas (hard spheres) of N (monatomic) particles at temperature T
#In a fixed-volume cube from -L:L, volume = 8L^3
#Can choose either solid or periodic walls


Ns = [10,10,10]           #number of atoms of different species [N1, N2, ...]
atom_types = ['He','Ar','Ne']


T = 300          #Kelvin, pick the desired temperature.
                #T <~ 1e7 gives viewable results
                #v_avg = sqrt(3kT/m)
                #Simulation does not accurately represent phase changes
                #T_boil (He-4) = 4.22K
                
Lnfac = 0.55     # >=0.55 (recommended), scales how big to make the box

g = 9.8          #m/s^2, =0 and =9.8 to see effects (minimal)

walls_type = 'periodic'     #'solid', or 'periodic'
n_clear_trail = 5000        #number of points after which to clear the trail
make_graphs = True
particle_collisions = True  #if false, they pass through, no way to Xfer energy
trail_color=color.red

the_rate = 3000
the_rate0 = the_rate
dt_min = 1e-15
k = 1.38064852e-23  #J/K, boltzmann constant

atom_dict = {'He':{'m':6.646e-27,'r':3.1e-10},'Ne':{'m':3.35e-26,'r':3.8e-10},
                'Ar':{'m':6.63e-26,'r':7.1e-10}} #mass in kg, radius in meters
rs = [atom_dict[a]['r'] for a in atom_types]
ms = [atom_dict[a]['m'] for a in atom_types]

#m, limit of cubic volume, volume = 8L^3: -L --> L all dimensions
#adjust L so it provides decent size
L = 0
for i, N in enumerate(Ns):
    L += Lnfac*rs[i]*N    #+= means L = L + rs[i]*N
    #1/(n-dimensions) = 1/3, with a factor of 2 and some comfort

canv_width=1000
canv_height=canv_width*(9/16)
scene=display(width=canv_width,height=canv_height,autoscale=False,range=1.8*L)

scene.caption = "Vary the animation rate: \n"
def set_rate(s):
    global the_rate
    the_rate = s.value

slider(min=the_rate0/10,value=the_rate0,max=the_rate0*10,length=0.5*canv_width,bind=set_rate)

scene.append_to_caption('\n\n')

#check if dt is too small
m_min = ms[0]
for m in ms:
    if m_min > m:
        m_min = m
dt = L/(100*sqrt(3*k*T/m_min))  #set dt: atoms don't move too far in one step
if dt < dt_min:
    Lfac = dt_min/dt    #greater than 1
    L *= Lfac
    dt = dt_min
    print('dt:',dt, 'L:',L, 'Lfac:',Lfac,'dt_orig:',dt_min/Lfac)
  

#draw cube: -L --> L
#periodic walls?
hf_thick=rs[0]
thick=hf_thick*2
the_color = color.white
if walls_type == 'periodic':
    opac = 0.25
    color_p = color.orange
    color_n = color.blue
else:
    opac = 1.0
    color_p = the_color
    color_n = the_color
L2 = 2*L
#the_walls = sphere(pos=vec(0,0,0),radius=1.5*L,color=color.green,opacity=0.25)

left_wall = box(pos=vec(-L-hf_thick,0,0),axis=vec(1,0,0),length=thick,
            width=L2,height=L2+2*thick,up=vec(0,1,0),color=color_n,opacity=opac)
right_wall = box(pos=vec(L+hf_thick,0,0),axis=vec(-1,0,0),length=thick,
            width=L2,height=L2+2*thick,up=vec(0,1,0),color=color_p,opacity=opac)
bottom_wall = box(pos=vec(0,-L-hf_thick,0),axis=vec(0,1,0),length=thick,
            width=L2,height=L2,color=color_n,opacity=opac)
top_wall = box(pos=vec(0,L+hf_thick,0),axis=vec(0,-1,0),length=thick,
            width=L2,height=L2,color=color_p,opacity=opac)
back_wall = box(pos=vec(0,0,-L-hf_thick),axis=vec(0,0,1),length=thick,
            width=L2+2*thick,height=L2+2*thick,color=color_n,opacity=opac)
#no front wall

def place_atoms(Ns, rs, L, method):
    #initialize position of objects: given L can fit all the atoms 
    colors = [color.green, color.orange, color.magenta, color.yellow, 
        color.cyan, vec(1,0,1)]
    
    atoms = []
    sfac = 1         #> 1 and < 1.5, shift factor
    if method == 'axes' or method == 'axis':
        #spread them out along one axis (x)
        ogn=vec(0,0,0)  #starting location, shift from here
        Ntot = 0
        shift = 0.25*rs[0]
        n = 0       #atom count
        Ntot = 0
        for i_N, N in enumerate(Ns):
            Ntot += N
            shift += sfac*rs[i_N]
            i=0
            while n < Ntot:
                if n%6 == 0:
                    shift += sfac*rs[i_N]
                if n%2 == 0:
                    sign=1
                else:
                    sign=-1

                if int(n/2)%3 == 0:
                    shift_vec=vec(shift,0,0)
                elif int(n/2)%3 == 1:
                    shift_vec=vec(0,shift,0)
                else:
                    shift_vec=vec(0,0,shift)

                if n == 0:
                    mt = True
                    if len(Ns) > 1:
                        the_color = colors[i_N%len(Ns)]
                    else:
                        the_color = trail_color
                else:
                    mt = False
                    the_color = colors[i_N%len(Ns)]
                #never run out of colors
                atoms.append(sphere(pos=ogn+sign*shift_vec,
                    radius=rs[i_N],color=the_color,make_trail=mt))
                if n == 0 and len(Ns) > 1:
                    atoms[-1].trail_color = trail_color
                n += 1  #increase running atom count by 1
                i += 1
                if n%6 == 0:
                    shift += sfac*rs[i_N]
                
    return atoms
    
    
def kickstart_atoms(atoms,Ns,ms,T):
    #give atoms momentum: mass and velocity
    k=1.38064852e-23    #J/K
    n = 0   #total atom count
    Ntot = 0
    for i_N, N in enumerate(Ns):
        Ntot += N
        while n < Ntot:
            atoms[n].m = ms[i_N]
            atoms[n].v = vector.random()
            n += 1
    
    #n will be the total number of atoms
    #Normalize by kinetic energy
    Ek=0        #joules
    for a in atoms:
        Ek += 0.5*a.m*mag2(a.v)
    #Need Ek=Ntot*1.5kT: scale all velocities to the average
    vfac = sqrt(1.5*n*k*T/Ek)
    #print('Ek:',Ek,'vfac:',vfac,'n:',n)
    #Ek = 0
    for a in atoms:
        a.v = a.v*vfac
        #print('vel:',a.v)
        #Ek += 0.5*a.m*mag2(a.v)    #to double check vfac is correct
    #print('Ek:',Ek, '= 1.5NkT =',1.5*n*k*T)    #passed
    
    return atoms

def check_walls(atoms,L,method):
    dp = 0
    if method == 'solid':
        for a in atoms:
            #left
            neg_lim = -L+a.radius
            pos_lim = L-a.radius
            neg_diff = a.pos - neg_lim*vec(1,1,1)
            pos_diff = a.pos - pos_lim*vec(1,1,1)
            #print('neg_diff:',neg_diff, 'pos_diff:',pos_diff)
            if neg_diff.x <= 0:      #left
                #shift back inside by too-far amount
                a.pos.x = neg_lim - neg_diff.x
                a.v.x = -a.v.x      #elastic bounce means v in normal dir flips
                dp += 2*a.m*abs(a.v.x)
            elif pos_diff.x >= 0:    #right
                #shift back inside by too-far amount
                a.pos.x = pos_lim - pos_diff.x
                a.v.x = -a.v.x      #elastic bounce means v in normal dir flips
                dp += 2*a.m*abs(a.v.x)
            if neg_diff.y <= 0:      #bottom
                #shift back inside by too-far amount
                a.pos.y = neg_lim - neg_diff.y
                a.v.y = -a.v.y      #elastic bounce means v in normal dir flips
                dp += 2*a.m*abs(a.v.y)
            elif pos_diff.y >= 0:    #top
                #shift back inside by too-far amount
                a.pos.y = pos_lim - pos_diff.y
                a.v.y = -a.v.y      #elastic bounce means v in normal dir flips
                dp += 2*a.m*abs(a.v.y)
            if neg_diff.z <= 0:      #back
                #shift back inside by too-far amount
                a.pos.z = neg_lim - neg_diff.z
                a.v.z = -a.v.z      #elastic bounce means v in normal dir flips
                dp += 2*a.m*abs(a.v.z)
            elif pos_diff.z >= 0:    #front
                #shift back inside by too-far amount
                a.pos.z = pos_lim - pos_diff.z
                a.v.z = -a.v.z      #elastic bounce means v in normal dir flips
                dp += 2*a.m*abs(a.v.z)
    elif method == 'periodic':
        #simplest way is to just move it, right at +/- L
        #other way is to create a particle while it straddles the wall

        for i, a in enumerate(atoms):
            neg_diff = a.pos - (-L*vec(1,1,1))
            pos_diff = a.pos - (L*vec(1,1,1))
            if i == 0:
                if neg_diff.x<0 or neg_diff.y<0 or neg_diff.z<0 or \
                    pos_diff.x>0 or pos_diff.y>0 or pos_diff.z>0:
                        a.make_trail = False
            if neg_diff.x < 0:
                a.pos.x = L + neg_diff.x
            elif pos_diff.x > 0:
                a.pos.x = -L + pos_diff.x
            if neg_diff.y < 0:
                a.pos.y = L + neg_diff.y
            elif pos_diff.y > 0:
                a.pos.y = -L + pos_diff.y
            if neg_diff.z < 0:
                a.pos.z = L + neg_diff.z
            elif pos_diff.z > 0:
                a.pos.z = -L + pos_diff.z
            if i == 0:
                a.make_trail = True

    return [atoms,dp]
    
def check_particles(atoms):
    for i,a in enumerate(atoms[:-1]):
        for b in atoms[i+1:]:
            ab = a.pos - b.pos
            RaRb2 = (a.radius+b.radius)**2
            if mag2(ab) <= RaRb2:
                M = a.m + b.m
                #back up to where they should have contacted
                
                Mr2 = M*mag2(ab) #Mr2 = M*mag2(ab)
                rv = dot(ab,b.v - a.v)
                vba = mag(b.v - a.v)
                rv0 = rv/vba
                b2 = mag2(ab) - rv0**2
                del_t = (sqrt(RaRb2 - b2) - rv0)/vba
                a.pos = a.pos - a.v*del_t
                b.pos = b.pos - b.v*del_t
                #compute new velocities, acquired through cm calculation
                a.v = a.v + 2*(b.m/Mr2)*rv*ab       #verified
                b.v = b.v - 2*(a.m/Mr2)*rv*ab
                
                
                #how far out from impact they should be
                a.pos = a.pos + a.v*del_t    #scoot them out
                b.pos = b.pos + b.v*del_t    #from impact

    return atoms
            

#initialize position of objects
atoms=place_atoms(Ns,rs,L,'axes')

#initialize velocity of objects
atoms=kickstart_atoms(atoms,Ns,ms,T)

#To see effect of gravity, perhaps on center of mass in y-dir
for a in atoms:
    a.a = vec(0,-g,0)

amag = mag(a)
scene.waitfor('click')

#loop
t = 0
n = 0

#graphs
if make_graphs:
    n_plot=50
    Ekscale = 0.5*(atoms[-1].m*mag2(atoms[-1].v) + atoms[-2].m*mag2(atoms[-2].v))
    #Ekscale < 1
    Ekscale = int(log(Ekscale)/log(10))

    #Ntot=sum(Ns)
    #P_ideal = Ntot*k*T/(8*L**3)
    cw2 = 0.45*canv_width
    ch2 = 0.45*canv_height
    cm_g = graph(width=cw2,height=ch2,align="left",ymin=-1,ymax=1)
    v_g = graph(width=cw2, height=ch2,align="right")
    E_g = graph(width=cw2,height=ch2,align="left")
    p_g = graph(width=cw2,height=ch2,align="right")
    #P_g = graph(width=cw2,height=ch2,align="left")
    cmx_s = series(graph=cm_g,label="Center-of-Mass-X/L", color=color.red)
    cmy_s = series(graph=cm_g,label="Center-of-Mass-Y/L", color=color.green)
    cmz_s = series(graph=cm_g,label="Center-of-Mass-Z/L", color=color.blue)    
    v_s = series(graph=v_g, label="tracer atom's speed (m/s)", color=trail_color)
    E_s = series(graph=E_g,label='Ek_tot (10^'+str(Ekscale)+' J)',color=color.black)
    px_s = series(graph=p_g,label='ptot_x (kg m/s)',color=color.red)
    py_s = series(graph=p_g,label='ptot_y (kg m/s)',color=color.green)
    pz_s = series(graph=p_g,label='ptot_z (kg m/s)',color=color.blue)
    #P_s = series(graph=P_g,label='Pressure(Pa)=?'+str(P_ideal),color=color.red)#avg for n=n_plot

    Ekscale = 10**-Ekscale
#P = 0
while t < 1:
    rate(the_rate)
    if make_graphs:
        if n%n_plot==0:
        #    P *= ((n_plot-1)*dt/(24*L**2))   #average pressure during dt*n_plot
            #plot energy and momentum

            Ektot = 0
            ptot = vec(0,0,0)
            #ptot_mag = 0
            cm = vec(0,0,0)
            Mtot = 0       
            for i, N in enumerate(Ns):
                Mtot += ms[i]*N
            for a in atoms:
                Ektot += 0.5*a.m*mag2(a.v)
                ptot = ptot + a.m*a.v
                cm = cm + a.m*a.pos
            cm = cm/(L*Mtot)
            cmx_s.plot(t,cm.x)
            cmy_s.plot(t,cm.y)
            cmz_s.plot(t,cm.z)
            v_s.plot(t,mag(atoms[0].v))
            E_s.plot(t,Ektot*Ekscale)
            px_s.plot(t,ptot.x)
            py_s.plot(t,ptot.y)
            pz_s.plot(t,ptot.z)
        #    P_s.plot(t,P)
        #    P = 0
        #else:
        #    P += dp
    
    #plot averagepressure on walls
    #plot speed distribution

    #make them move
    for i,a in enumerate(atoms):
        a.pos = a.pos + a.v*dt
        if amag != 0:
            a.v = a.v + a.a*dt
    
    #check for interparticle collisions
    if particle_collisions:
        atoms = check_particles(atoms)
    
    #check for wall collisions
    [atoms,dp] = check_walls(atoms,L,walls_type)

    if n%n_clear_trail == 0:
        atoms[0].clear_trail()
        
    n += 1
    t = n*dt
