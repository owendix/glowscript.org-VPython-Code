GlowScript 2.2 VPython
#Shows the image and two principal rays for an object for different 
#thin geometric optics: convex lens, concave lens, convex mirror, concave mirror
h_o_fac = 0.666     #0.666, good visibility, smaller for better thin optic approx

show_oi = True

#to invert the colors on a Mac: Control-Option-Command-8

f_mag = 1       #magnitude of focus
f = f_mag
the_range = 1.95*f_mag
ax_r = 0.01*the_range
vis_range = 3*the_range
optic_types = ['convex lens','concave lens','convex mirror','concave mirror']
optic_type = optic_types[0]

optic_blue = vec(0.75,0.75,1)
obj_color = color.cyan
img_color = color.orange
optic_opac = 0.5
image_opac = 0.85
vray_opac = 0.333
silver = vec(0.75,0.75,0.75)
f_color=color.red
c_color=color.blue
ray_color=color.white
is_white = ray_color.equals(color.white)
ray_r = 0.5*ax_r

#create scene
cw = 1200
ch = (1/3)*cw
#locked scene
scene=display(width=cw, height=ch, userzoom=False,userspin=False,autoscale=False)
scene.range=the_range
scene.title='Drag object to see image location for different types of optics\n\n'

#create center axis

ax = cylinder(pos=vec(-vis_range,0,0),axis=vec(2*vis_range,0,0),
            radius=ax_r,color=color.white, opacity=0.5)
vax = cylinder(pos=vec(0,-vis_range,0,0),axis=vec(0,2*vis_range,0),
            radius=0.5*ax_r,color=color.white, opacity=0.5)
optic = 0
object = 0
image = 0
rays=[[],[],[]]  #rays = [[ray1], [ray2], [ray3]]
fs = []
c_curve = 0
optic_tails = []
o_label=0
i_label=0
o = 2*f_mag    #must be positive (left side), set initial object location
i = o*f_mag/(o-f_mag)
ang0 = 0
op_h = 3.5*f_mag

def make_optic():
    global optic, f, optic_tails, ang0, fs, c_curve
    #center optic at origin

    if len(fs) > 0:
        for q in fs:
            q.visible=False
        fs=[]
        optic.visible = False
        for q in optic_tails:
            q.visible=False
        optic_tails = []
    
    if c_curve != 0:
        c_curve.visible = False

    if 'lens' in optic_type:
        for pn in [-1,1]:
            fs.append(sphere(pos=vec(pn*f_mag,0,0),radius=1.5*ax_r,color=f_color))
        optic = cylinder(pos=vec(0,-0.5*op_h,0), axis=vec(0,op_h,0),
                    radius=ax_r, color=optic_blue, opacity=optic_opac)
        l_thing = 10*ax_r
        if 'vex' in optic_type: #convex/converging
            #draw top and bottom things signifying
            for fp in fs:
                optic_tails.append(cylinder(pos=optic.pos,
                    axis=l_thing*hat(fp.pos-optic.pos),radius=ax_r,
                    color=optic_blue,opacity=optic_opac))
                optic_tails.append(cylinder(pos=optic.pos+optic.axis,
                    axis=l_thing*hat(fp.pos-(optic.pos+optic.axis)),radius=ax_r,
                    color=optic_blue,opacity=optic_opac))
            f = f_mag
        else:       #concave/diverging lens
            #draw top and bottom things
            for fp in fs:
                optic_tails.append(cylinder(pos=optic.pos,
                    axis=-l_thing*hat(fp.pos-optic.pos),radius=ax_r,
                    color=optic_blue, opacity=optic_opac))
                optic_tails.append(cylinder(pos=optic.pos+optic.axis,
                    axis=-l_thing*hat(fp.pos-(optic.pos+optic.axis)),radius=ax_r,
                    color=optic_blue, opacity=optic_opac))
            f = -f_mag
    else:
        crvlist = []
        if 'vex' in optic_type:
            f = -f_mag
        else:
            f = f_mag
        
        fs.append(sphere(pos=vec(-f,0,0),radius=1.5*ax_r,color=f_color))
        
        c_ogn = vec(-2*f,0,0)
        c_curve = sphere(pos=c_ogn, radius=1.5*ax_r, color=c_color,visible=True)
        npts = 15       #must be odd for symmetry
        c_r = -c_ogn.x  #positive when ogn on left, neg when on right
        ang = asin(op_h/(2*abs(c_r)))   #will be positive
        ang0 = ang
        dang = -2*ang/npts  #will be negative
        while ang > -ang0:
            crvlist.append(c_ogn + c_r*vec(cos(ang),sin(ang),0))
            ang += dang
        optic=curve(pos=crvlist,radius=ax_r,color=silver)

def make_object():
    global object, o
    #object must be on negative side
    h_o = h_o_fac*f_mag
    if o < 0:
        o = -o
    elif o > 4*f_mag:
        o = 4*f_mag
    object = arrow(pos=vec(-o,0,0),axis=vec(0,h_o,0),shaftwidth=0.1*abs(h_o),
                    color=obj_color)
                    


def make_image(is_default):
    global image, rays, o, i, o_label, i_label
    o = -object.pos.x
    i = o*f/(o-f)    
    h_i = -(i/o)*object.axis.y
    h_iy = -(i/o)*object.pos.y
    if 'lens' in optic_type:
        try:
            image.pos = vec(i,h_iy,0)
            image.axis = vec(0,h_i,0)
            sw = 0.1*abs(h_i)
            image.shaftwidth = sw
            image.headwidth = 2*sw
            image.headlength = 3*sw
        except:
            sw = 0.1*abs(h_i)
            image = arrow(pos=vec(i,h_iy,0),axis=vec(0,h_i,0),
                    shaftwidth=sw, headwidth = 2*sw, headlength=3*sw,
                    color=img_color, opacity=image_opac)
    else:
        try:    #if image != 0
            image.pos=vec(-i,h_iy,0)
            image.axis=vec(0,h_i,0)
            sw = 0.1*abs(h_i)
            image.shaftwidth = sw
            image.headwidth = 2*sw
            image.headlength = 3*sw
        except:
            sw = 0.1*abs(h_i)
            image = arrow(pos=vec(-i,h_iy,0),axis=vec(0,h_i,0),
                    shaftwidth=sw, headwidth = 2*sw, headlength=3*sw,
                    color=color.orange, opacity=image_opac)
    if show_oi:
        if not is_default:
            o_label.text='o='+'{0:.2f}'.format(o)+'f'
            i_label.text='i='+'{0:.2f}'.format(i)+'f'
        else:
            o_label=label(text='o='+'{0:.2f}'.format(o)+'f',pos=vec(-5.75*f_mag,-0.45*op_h,0),
                    align='left',height=16, font='monospace',
                    box=False,line=False,xoffset=0,yoffset=0,opacity=0,
                    color=obj_color,visible=True)
            i_label=label(text='i='+'{0:.2f}'.format(i)+'f',pos=vec(-5.75*f_mag,-0.5*op_h,0),
                    align='left',height=16, font='monospace',
                    box=False,line=False,xoffset=0,yoffset=0,opacity=0,
                    color=img_color,visible=True)
        

make_optic()
make_object()
make_image(True)
            

drag = False

def inside_object(a_pos):
    #approximate arrow as a rectangle, vertically aligned
    #position components of [right, top, left, bottom]
    wd2 = 0.5*object.shaftwidth
    rect=[object.pos.x+wd2, 
        object.pos.y+object.axis.y,
        object.pos.x-wd2,
        object.pos.y]
    #print(wd2, rect, a_pos)
    if a_pos.x >= rect[0] or a_pos.x <= rect[2] or \
            a_pos.y >= rect[1] or a_pos.y <= rect[3]:
        return False
        
    return True
    
diff_vec0 = vec(0,0,0)

scene.bind("mousedown", def ():
    global drag, diff_vec0
    drag = inside_object(scene.mouse.pos)
    diff_vec0 = object.pos - scene.mouse.pos
    #print(drag)
)

scene.bind("mousemove", def ():
    global drag, object
    eps=0.1*abs(object.shaftwidth)
    if drag and object.pos.x <= -eps: # mouse button is down
        new_pos = scene.mouse.pos + diff_vec0
        if new_pos.x < -eps:
            object.pos = new_pos
        if object.pos.x > -eps:
            object.pos = object.pos - vec(eps,0,0)
            drag = False
        make_image(False)
        move_rays()
)

scene.bind("mouseup", def ():
    nonlocal drag
    drag = False
)


def set_optic(optic_radio):
    global optic_type, convexlens_radio, concavelens_radio, convexmirror_radio,\
        concavemirror_radio
    the_radios = [convexlens_radio,concavelens_radio,convexmirror_radio,
        concavemirror_radio]
    #radio seems to automatically toggle, next line toggles it back
    optic_radio.checked = not optic_radio.checked
    if 'lens' in optic_radio.text:
        if 'vex' in optic_radio.text:   #convex lens
            if convexlens_radio.checked:
                return  #it was already checked, do nothing
            for r in the_radios:
                r.checked = False
            optic_radio.checked = True
            optic_type=optic_types[0]               
            
        else:                           #concave lens
            if concavelens_radio.checked:
                return
            for r in the_radios:
                r.checked = False
            optic_radio.checked = True
            optic_type=optic_types[1]                

    else:       #mirror
        if 'vex' in optic_radio.text:   #convex mirror
            if convexmirror_radio.checked:
                return
            for r in the_radios:
                r.checked = False
            optic_radio.checked = True
            optic_type=optic_types[2]                

        else:                           #concave mirror
            if concavemirror_radio.checked:
                return
            for r in the_radios:
                r.checked = False
            optic_radio.checked = True
            optic_type=optic_types[3]
            
    make_optic()
    make_image()
    if show_rays[0]:
        make_ray1(bray1,from_default=True)
    if show_rays[1]:
        make_ray2(bray2,from_default=True)
            
        
#4 radios
scene.caption='\n<b>Select optic type: <b\>\n\n'
convexlens_radio = radio(text='convex lens',checked='vex l' in optic_type,
    bind=set_optic)
scene.append_to_caption('\t\t')
concavelens_radio = radio(text='concave lens',checked='cave l' in optic_type,
    bind=set_optic)
scene.append_to_caption('\t\t')
convexmirror_radio = radio(text='convex mirror',
    checked='vex m' in optic_type,bind=set_optic)
scene.append_to_caption('\t')
concavemirror_radio = radio(text='concave mirror',
    checked='cave m' in optic_type,bind=set_optic)


#PRINCIPAL RAYS

show_rays = [False, False]   #always start false, necessary for my code

scene.append_to_caption('\n\n<b>Draw Principal Rays: <b\>\n\n')

def make_ray1(b,from_default=False):
    global rays, show_rays
    #along axis, through/away from focus
    if not from_default:
        show_rays[0] = not show_rays[0]
        if 'Show' in b.text:
            b.text = ' Hide Ray 1 '
        else:
            b.text = 'Show Ray 1'
    if show_rays[0]:
        for q in rays[0]:
            q.visible = False
        rays[0] = []   #will be true to start
        tip=object.pos + object.axis
        rays[0].append(cylinder(pos=tip,
                axis=vec(-tip.x,0,0),radius=ray_r,color=ray_color))
        optic_impact = vec(0, rays[0][0].pos.y, 0)
        if 'lens' in optic_type:
            if 'vex' in optic_type:
                #fs[1] is always on the right side of the screen for lenses
                dir_real = vis_range*(fs[1].pos - optic_impact)
            else:
                dir_real = -vis_range*(fs[0].pos - optic_impact)
        else:           #mirror
            dir_real = vis_range*(fs[0].pos - optic_impact)
            if 'vex' in optic_type:
                dir_real = -dir_real
        rays[0].append(cylinder(pos=optic_impact,axis=dir_real,
                        radius=ray_r,color=ray_color))
        if i < 0:   #virtual image
            rays[0].append(cylinder(pos=optic_impact,axis=-dir_real,
                    radius=ray_r,color=ray_color,opacity=vray_opac))
    else:
        for q in rays[0]:
            q.visible = False
        rays[0] = []
    
#button ray1
bray1 = button(text='Show Ray 1', bind=make_ray1)
scene.append_to_caption('\t\t')

def make_ray2(b, from_default=False):
    global rays, show_rays
    #at optic center, undeflected or equal angles
    if not from_default:
        show_rays[1] = not show_rays[1]
        if 'Show' in b.text:
            b.text = ' Hide Ray 2 '
        else:
            b.text = 'Show Ray 2'
    if show_rays[1]:
        for q in rays[1]:
            q.visible = False
        rays[1] = []   #will be true to start
        tip=object.pos + object.axis
        rays[1].append(cylinder(pos=tip,
                axis=-tip,radius=ray_r,color=ray_color))
        optic_impact = vec(0, 0, 0)
        if 'lens' in optic_type:
            dir_real = -vis_range*hat(tip)
        else:           #mirror
            dir_real = vis_range*hat(vec(tip.x,-tip.y,tip.z) - optic_impact)
        rays[1].append(cylinder(pos=optic_impact,axis=dir_real,
                        radius=ray_r,color=ray_color))
        if i < 0:   #virtual image
            rays[1].append(cylinder(pos=optic_impact,axis=-dir_real,
                    radius=ray_r,color=ray_color,opacity=vray_opac))
    else:
        for q in rays[1]:
            q.visible = False
        rays[1] = []
    
#button ray2
bray2 = button(text='Show Ray 2', bind=make_ray2)
scene.append_to_caption('\t\t')


thin_optic=label(text='thin optics approximation: rays are drawn to vertical center'
            +' line and optics are assumed sufficiently tall',
            pos=vec(0,-0.5*op_h,0),align='center',color=color.white,box=False,
            line=False,height=14,opacity=0,yoffset=-1)

#move principal rays
def move_rays():
    global rays
    if show_rays[0]:
        tip=object.pos + object.axis
        rays[0][0].pos=tip
        rays[0][0].axis=vec(-tip.x,0,0)
        optic_impact = vec(0, rays[0][0].pos.y, 0)
        if 'lens' in optic_type:
            if 'vex' in optic_type:
                #fs[1] is always on the right side of the screen for lenses
                dir_real = vis_range*hat(fs[1].pos - optic_impact)
            else:
                dir_real = -vis_range*hat(fs[0].pos - optic_impact)
        else:           #mirror
            dir_real = vis_range*(fs[0].pos - optic_impact)
            if 'vex' in optic_type:
                dir_real = -dir_real
        rays[0][1].pos=optic_impact
        rays[0][1].axis=dir_real
        if i < 0:   #virtual image
            if len(rays[0]) == 2:
                
                rays[0].append(cylinder(pos=optic_impact,axis=-dir_real,
                    radius=ray_r,color=ray_color,opacity=vray_opac))
            else:
                rays[0][2].pos=optic_impact
                rays[0][2].axis=-dir_real
        elif i > 0:
            if len(rays[0]) > 2:
                rays[0][2].visible=False
                rays[0] = rays[0][:-1]
                
    if show_rays[1]:
        tip=object.pos + object.axis
        rays[1][0].pos=tip
        rays[1][0].axis=-tip
        if 'lens' in optic_type:
            dir_real = -vis_range*hat(tip)
        else:           #mirror
            dir_real = vis_range*hat(vec(tip.x,-tip.y,tip.z))
        rays[1][1].pos=vec(0,0,0)
        rays[1][1].axis=dir_real
        if i < 0:   #virtual image
            if len(rays[1]) == 2:
                rays[1].append(cylinder(pos=vec(0,0,0),axis=-dir_real,
                    radius=ray_r,color=ray_color,opacity=vray_opac))
            else:
                rays[1][2].pos=vec(0,0,0)
                rays[1][2].axis=-dir_real
        elif i > 0:
            if len(rays[1]) > 2:
                rays[1][2].visible=False
                rays[1] = rays[1][:-1]


def wavelength_to_rgb(wavelength):
    
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
    lam_r = 750
    lam_b = 380
    gamma=0.8
    
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
    
    return R, G, B

scene.append_to_caption('\n\n<b>Vary the principal ray color (wavelength): <b\>\n\n')

def vary_color(lc):
    global is_white, ray_color, rays, lc_radio
    
    if 'white' in lc.text:  #its the radio
        is_white = lc_radio.checked
        #lc_radio.checked = is_white
    else: #its the slider
        #print('its the slider')    #works
        is_white = False
        lc_radio.checked = False
    
    #change the light_color variable
    if not is_white:
        lc_value = lc_slider.value
        R,G,B = wavelength_to_rgb(lc_value*1e-9)
        ray_color = vec(R,G,B)

    else:
        ray_color = color.white
    #change the colors of the light
    for q in rays:
        for r in q:
            r.color=ray_color
            
#wavelength in nanometers
lc_slider=slider(value=585,min=380,max=750,step=1,length=0.4*cw,bind=vary_color)
lc_slider.text='color'  #append for checking purpose

scene.append_to_caption('\t')

lc_radio = radio(text='white',checked=is_white, bind=vary_color)

