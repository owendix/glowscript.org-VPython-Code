GlowScript 2.1 VPython
"""
Simulates Laue X-ray back-scattering diffraction pattern for CeIrIn5,
along both basis vectors
Only simulates diffraction spikes in the half solid 
angle opposite the incident beam direction (Laue back scattering)
"""
#default, first one will be flipped (easiest way)
which_s0 = 1        #0-->a-axis, 1 --> c-axis (leave this as a default here)
stereo_proj = False     #stereographically projects (like a real x-ray), sphere if false
false_color = False      #maps x-ray wavelength into visible, white if false

#constants
h = 6.62607004e-34  #Js, for the energy of a photon
v = 299792458       #m/s, speed of light (X-rays)
e = 1.60217662e-19  #C, charge of electron
dV = 3.5e4              #Volts: voltage difference for electron producing X-rays
lam_min = h*v/(e*dV)    #smallest wavelength limited by max energy of electrons
#lam_min = 0.3542e-10   #m, this is the result of previous values
lam_max = 1.38059e-10   #m, not sure what determines this
N = 2500              #number of diffracted beam spots to print
params = [0 for q in range(N)]
beams = [0 for q in range(N)]

def wavelength_to_rgb(wavelength,lam1=lam_min,lam2=lam_max):
    
    '''This converts a given wavelength of light to an 
    approximate RGB color value, with a false color spectrum
    This scale the wavelength to the visible light spectrum
    
    Based on code by http://www.noah.org/wiki/Wavelength_to_RGB_in_Python,
    which was based on code by Dan Bruton
    http://www.physics.sfasu.edu/astro/color/spectra.html
    '''
    #Scale it so it doesn't have to be in the visible spectrum (false color)
    
    #Convert meters to nanometers
    wavelength = wavelength*1e9
    lam1 = lam1*1e9
    lam2 = lam2*1e9
    lam_r = 750
    lam_b = 380
    gamma=0.8
    if wavelength > lam_r or wavelength < lam_b:
        wavelength = lam_b+(lam_r-lam_b)*((wavelength-lam1)/(lam2-lam1))
    
    if wavelength >= 380 and wavelength <= 440:
        attenuation = 0.3 + 0.7 * (wavelength - 380) / (440 - 380)
        R = ((-(wavelength - 440) / (440 - 380)) * attenuation) ** gamma
        G = 0.0
        B = (1.0 * attenuation) ** gamma
    elif wavelength >= 440 and wavelength <= 490:
        R = 0.0
        G = ((wavelength - 440) / (490 - 440)) ** gamma
        B = 1.0
    elif wavelength >= 490 and wavelength <= 510:
        R = 0.0
        G = 1.0
        B = (-(wavelength - 510) / (510 - 490)) ** gamma
    elif wavelength >= 510 and wavelength <= 580:
        R = ((wavelength - 510) / (580 - 510)) ** gamma
        G = 1.0
        B = 0.0
    elif wavelength >= 580 and wavelength <= 645:
        R = 1.0
        G = (-(wavelength - 645) / (645 - 580)) ** gamma
        B = 0.0
    elif wavelength >= 645 and wavelength <= 750:
        attenuation = 0.3 + 0.7 * (750 - wavelength) / (750 - 645)
        R = (1.0 * attenuation) ** gamma
        G = 0.0
        B = 0.0
    else:
        R = 0.0
        G = 0.0
        B = 0.0
    
    return (R, G, B)


def set_dir(ds0):
    if ds0 == 0:
        s0 = vec(-1,0,0)    #only these two directions, future script depends on it
        the_axis = 'a'
    else:
        s0 = vec(0,0,-1)
        the_axis = 'c'
        
    hkl_lim = 20
    return [s0, the_axis, hkl_lim]

[s0, the_axis, hkl_lim] = set_dir(which_s0)
    
canv_w = 900
canv_h = (9/16)*canv_w
scene=display(width=canv_w,height=canv_h,up=vec(0,1,0),center=vec(0,0,0))
#scene.camera.pos = -D*s0


def project_Laue(proj_radio,from_axis=False):
    global beams, stereo_proj, false_color
    
    lenFact = 225
    D = 3.4772e-2       #m, distance between screen and sample, at (0,0,0)
    L = lenFact*D #length of cylinders, D is distance to camera
    Iscale = 10000        #print and experiment to find a reasonable number
    if not from_axis:
        if 'axis' not in proj_radio:
            #Changing the color 
            if 'Color' in proj_radio.text or 'White' in proj_radio.text:
                #print('In Color part')
                #its the false_color change, don't regenerate beams
                #just change the colors
                if false_color and 'White' in proj_radio.text:
                    false_color=False
                    proj_radio_color.checked=False
                    proj_radio_white.checked=True
                elif (not false_color) and ('Color' in proj_radio.text):
                    false_color=True
                    proj_radio_color.checked=True
                    proj_radio_white.checked=False
                else:       #same one clicked, don't regenerate
                    return
            else:   #called by the stereographic or spherical projection
                if stereo_proj and 'Stereo' in proj_radio.text:  #Same projection
                    return
                elif (not stereo_proj) and 'Spher' in proj_radio.text: #same proj
                    return
                else:
                    if stereo_proj:   #change to a spherical proj
                        stereo_proj = False
                        proj_radio_sphere.checked=True
                        proj_radio_stereo.checked=False
                        scene.ambient=color.gray(0.1)
                    elif not stereo_proj:   #change to a stereo proj
                        stereo_proj = True
                        proj_radio_sphere.checked=False
                        proj_radio_stereo.checked=True
                        scene.ambient=color.gray(0.9)

    #print(len(params))
    #N: number of beam spots to print
    rangeFact = 5 
    radDenom = 50
    if stereo_proj:
        scene.range=rangeFact*L
    else:
        scene.range=1.2*L
    beam_rad_max = L/radDenom
    for b in beams:
        try:
            b.visible=False
        except:
            pass
            
    i=0
    for p in params:
        #opac = p[2]/I_max
        opac = 1 - 1/(1+p[2]/Iscale)
        if opac > 0.2:
            #print(p[0])
            #if opac > 1:
            #    opac = 1
            #opac -= 0.5     #shift
            #if opac < 0:
            #    opac = 0
            beam_rad = beam_rad_max*opac
            if false_color:
                R,G,B = wavelength_to_rgb(p[1])
                the_color=vec(R,G,B)
            else:
                the_color=color.white
            if not stereo_proj:
                try:
                    beams[i].pos=L*p[0]
                    beams[i].radius=beam_rad
                    beams[i].color=the_color
                    beams[i].visible=True
                except:
                    beams[i] = sphere(pos=L*p[0],radius=beam_rad,color=the_color,
                        visible=True)
                i += 1
            elif which_s0 == 0: #s0 = vec(-1,0,0)
                #p has unit radius, comes from s
                if abs(degrees(diff_angle(p[0],-s0))) < 85: #avoids division by 0
                    try:
                        beams[i].pos=L*vec(1,p[0].y/p[0].x, p[0].z/p[0].x)
                        beams[i].radius=beam_rad
                        beams[i].color=the_color
                        beams[i].visible=True
                    except:
                        beams[i] = sphere(pos=L*vec(1,p[0].y/p[0].x,p[0].z/p[0].x),
                                radius=beam_rad,color=the_color,visible=True)
                    i += 1
            else:               #s0 = vec(0,0,-1)    
                if abs(degrees(diff_angle(p[0],-s0))) < 85:   #avoids division by 0
                    try:
                        beams[i].pos=L*vec(p[0].x/p[0].z,p[0].y/p[0].z,1)
                        beams[i].radius=beam_rad
                        beams[i].color=the_color
                        beams[i].visible=True
                    except:
                        beams[i] = sphere(pos=L*vec(p[0].x/p[0].z,p[0].y/p[0].z,1),
                                radius=beam_rad,color=the_color,visible=True)
                    i += 1
    #end of generate_Laue

def generate_Laue(axis_radio,from_default=False):
    #generate the Laue diffraction pattern
    global which_s0, s0, params, params
    
    if not from_default:
        if which_s0 == 0 and ('c-axis' in axis_radio.text):      #prior to click
            which_s0 = 1
            dir_radio_c.checked=True
            dir_radio_a.checked=False
            scene.lights=[distant_light(direction=vec(0.22, 0.44, 0.88),color=color.gray(0.8)),
                distant_light(direction=vec(-0.88, -0.22, -0.44), color=color.gray(0.3))]
        elif which_s0 == 1 and ('a-axis' in axis_radio.text):
            which_s0 = 0
            dir_radio_c.checked=False
            dir_radio_a.checked=True
            scene.lights=[distant_light(direction=vec(0.88,0.44,0.22),color=color.gray(0.8)),
                distant_light(direction=vec(-0.44,-0.22,-0.88),color=color.gray(0.3))]
        else:
            return

    [s0, the_axis, hkl_lim] = set_dir(which_s0)
    scene.forward=s0
    
    
    #From Cullity, common limitation of X-ray photo paper: 
    #lam_max = 4e9          #m,
    #CeIrIn5 tetragonal crystal structure parameters
    a = 4.6662e-10          #m
    c = 7.5168e-10          #m
    #Iridium and Indium atoms at sites:
    ZCe = 58    # of protons/electrons, roughest approx for atomic form factor
    ZIr = 77    #assumes all electrons concentrated on the cell
    ZIn = 49    #I will ignore temperature effects: the Debye-Waller factor
    #rbasis: Ce, Indium, Iridium0, Iridium1, Iridium2, Iridium3, Iridium4
    p = 0.3035
    w = 0.5
    rbasis = [vec(0,0,0), vec(0,0,w), vec(w,w,0), vec(0,w,p), vec(w,0,p),
                vec(0,w,(1-p)), vec(w,0,(1-p))]
    Zbasis = [ZCe, ZIr, ZIn, ZIn, ZIn, ZIn, ZIn]
    #s0 is the incident beam direction
    #s will be the diffracted beam UNIT vector
    
    
    #crystal at origin: vec(0,0,0)
    #if falseColor:
    #    print('X-ray wavelength is rescaled into the visible spectrum')
    M = 100              #number of unit cells in intensity calculation
    #smaller hkl-absolute values are preferred
    if s0.x == -1:
        #Extra constraint based on incident beam direction
        dir_min = ceil(a/lam_max)
        dir_max = int(2*a/lam_min)
    elif s0.z == -1:
        #Extra constraint based on incident beam direction
        dir_min = ceil(c/lam_max)
        dir_max = int(2*c/lam_min)
    
    if hkl_lim < dir_max:
        dir_max = hkl_lim
    dir_list = [q for q in range(dir_min,dir_max+1)]
    #Do the non-weird ones
    ndir_list=[]
    for q in range(0,hkl_lim+1):
        for pn in [1,-1]:
            ndir_list.append(pn*q)
    ndir_list = ndir_list[1:]

    k_list = ndir_list.copy()

    #Starts with lowest abs-val of h,k,l's so its likeliest most intense
    """
    This takes too long in the multi-loop (tried it)
    def gcd(a,b):
        #efficient algorithm for finding greatest common divisor
        #I will apply this twice because gcd(a,b,c)=gcd(a,gcd(b,c))=gcd(gcd(a,b),c)
        while b != 0:
            t = b
            b = a%b
            a = t
        return a
    #May be faster to add the points than to check for gcd
    #find smallest val of h,k,l
    #q = gcd(h,gcd(k,l))
    #if q > 1:
    #    [h,k,l] = [h/q,k/q,l/q]
    """
    """
    #Conditions for a successful beam:
    1) mag2(s) == 1 : I double checked and this is satisfied
        --> s.x**2 + s.y**2 + s.z**2 == 1
    2) dot(s,s0) < 0
        if s0 = vec(-1,0,0) --> 1 > s.x > 0 (unit vec)
            --> 2a/lam > h > a/lam
                --> Use to narrow h: int(2a/lam_min) > h > int(a/lam_max)
                --> Then pick k, l and use to find lam exactly for that choice
                --> Then use as a better condition for h??? (always consistent?)
                --> Ah, its quick. Might as well do it.
        if s0 = vec(0,0,-1) --> 1 > s.z > 0 (unit vec)
            --> 2c/lam > l > c/lam
                --> Use to narrow l: int(2a/lam_min) > l > int(a/lam_max)
                --> Then pick h, k and use to find lam exactly for that choice
                --> Then use as a better condition for l
    H = vec(h/a, k/a, l/c)
    lam = -2*dot(H,s0)/mag2(H)
    s - s0 = lam*H = lam*vec(h/a, k/a, l/c)
    """
    n = 0
    for dir in dir_list:
        if n > N:
            break
        #set vecs equal
        if s0.x == -1:
            h = dir
        elif s0.z == -1:
            l = dir
        for k in k_list:
            if n > N:
                break
            for ndir in ndir_list:
                if s0.x == -1:
                    l = ndir
                elif s0.z == -1:
                    h = ndir
                H = vec(h/a, k/a, l/c)
                lam = -2*dot(H,s0)/mag2(H)
                if lam < lam_min or lam > lam_max:
                    continue
                s = lam*H + s0
                #make sure its just that pointed at the screen
                if dot(s,s0)>=0:
                    continue
                #need all of this, especially if I want to adjust brightness
                #Calculate intensity
                I_lattice = 1
                if h != 0:
                    I_lattice *= (sin(M*h/2)/sin(h/2))**2
                if k != 0:
                    I_lattice *= (sin(M*k/2)/sin(k/2))**2
                if l != 0:
                    I_lattice *= (sin(M*l/2)/sin(l/2))**2
                if I_lattice == 0:
                    continue
                struct_factor = vec(0,0,0)
                for i,b in enumerate(rbasis):
                    cb = cos(2*pi*dot(b,vec(h,k,l)))
                    sb = sin(2*pi*dot(b,vec(h,k,l)))
                    struct_factor = struct_factor + Zbasis[i]*vec(cb,sb,0)
                I = I_lattice*mag2(struct_factor)
                if I == 0:
                    continue
                n += 1
                if n > N:
                    break
                #print(h,k,l,'struct_factor:',struct_factor,'\n')
                #print(I_lattice,mag2(struct_factor))
                params[n-1]=[s,lam,I,vec(h,k,l)]
    
    #params acquired
    #Similar to direction, pick other
    if stereo_proj:
        project_Laue(proj_radio_stereo,from_axis=True)
    else:
        project_Laue(proj_radio_sphere,from_axis=True)
    #end of generate Laue

scene.caption='One of CeIrIn_5\'s a-axis is vertical\n\n'
scene.append_to_caption('Beam incident along:\t\t')
scene.append_to_caption('Projection Type:\t\t\t\t\t')
scene.append_to_caption('Beam color\n')
dir_radio_a = radio(bind=generate_Laue,checked=(which_s0 == 0),pos=scene.caption_anchor,text='<b>a-axis<b\>')
scene.append_to_caption('\t\t\t\t')
proj_radio_sphere = radio(bind=project_Laue, checked=(not stereo_proj),text='<b>Spherical Projection<b\>')
scene.append_to_caption('\t\t\t')
proj_radio_color = radio(bind=project_Laue, checked=false_color,text='<b>Color Mapped<b\>')
scene.append_to_caption('\n')
dir_radio_c = radio(bind=generate_Laue,checked=(which_s0 == 1),text='<b>c-axis<b\>') # text to right of button
scene.append_to_caption('\t\t\t\t')
proj_radio_stereo = radio(bind=project_Laue,checked=stereo_proj,text='<b>Stereographic Projection<b\>')
scene.append_to_caption('\t\t')
proj_radio_white = radio(bind=project_Laue,checked=(not false_color),text='<b>White<b\>')
scene.append_to_caption('\n')

generate_Laue(dir_radio_c,from_default=True)

