GlowScript 2.1 VPython
#Create a force diagram from given forces
#and demonstrate tip-to-tail, tail-to-tip method
#Scale will be in Newtons: force-space
#Glowscript limits this to 4 vectors w/ magnitudes more (>6) without 
#magnitudes (creating text exceeds vertex #)

#System
sys = 'Ball'

show_magnitudes = True     #If false, don't show magnitudes
#The type and EXTERNAL Object Causing Force on System (not the system itself)
#Keep brackets unchanged, copy and paste brackets and edit for more forces
#This uses the Python dictionary built-in data structure
#No comma after last in list
#types:
#g --> gravity
#c --> contact
#e --> electric
#m --> magnetic
#e.g. [{{'type': 'g'}, {'obj': 'Earth'}},
#       {{'type': 'c'}, {'obj': 'Earth'}},      #last comma okay but not needed     
#       ]
#List in same order of Force values given below, no zero-length vectors
Fnames = [{'type': 'g', 'obj': 'Earth'},
        {'type': 'c', 'obj': 'Rope'}
        ]
                                               

#Force values in Newtons,view: +x=right, +y=up, +z=toward camera
#Components vec(Fx, Fy, Fz), same order as Fname above, no zero-length vectors
Fvals = [vec(0,-750,0),
        vec(0,300,0)
        ]                                       

#Set scale for force diagram
#Find largest and smallest forces
shift_rate = 1000
shift_fac = 100
sys_color = color.red
net_color = color.white
the_colors = [color.green, color.cyan, color.yellow, color.magenta, color.orange, 
                color.blue]
fmin = mag(Fvals[0])
fmax = mag(Fvals[0])
for f in Fvals:
    fmag = mag(f)
    if fmag < fmin:
        fmin = fmag
    elif fmag > fmax:
        fmax = fmag

rscale1 = 0.1
rscale2 = 0.5
if fmin < (rscale2*rscale1)*fmax:
    r_sys = rscale2*fmin
else:
    r_sys = rscale1*fmax
    
shaft_w = rscale2*r_sys
head_wl = 2*shaft_w
sys_name_scale = 1.5*fmax

canv_w = 1000
canv_h = (7/16)*canv_w
scene = display(width=canv_w, height=canv_h)
scene.title='CLICK TO ADVANCE'


def get_shift_hat(vec_orig, is_label):
    #How to shift certain things to avoid overlap
    dz = vec(0,0,1)  #shift outward for visibility
    if not is_label:
        return dz
    else:
        v = hat(cross(vec_orig, vec(0,0,-1)))
        #shift labels of more vertical forces less
        comp_y = dot(hat(vec_orig), vec(0,1,0))
        dy = vec(0,-0.5,0)*(1 + comp_y)
        scale = 0.5*(3 - 0.75*abs(comp_y))
        if vec_orig.x >= 0:
            if vec_orig.y >= 0:
                return scale*(-v + dz + dy)
            else:
                return scale*(v + dz + dy)
        else:
            if vec_orig.y >= 0:
                return scale*(-v + dz + dy)
            else:
                return scale*(v + dz + dy)

def shift_rects(r1, r2):
    #Determine if 2 rectangles overlap and return the best 
    #vertical shift of r1 to undo it
    #2 rectangles each in the following format:
    #[right,top,left,bottom]
    #rectangles known to be flatter: (top-bottom) < (right-left)
    
    shift_shortest = False
    fac = 1.1   #shift 10% farther
    #if false: always shift vertically away from zero (origin of f-diag)
    #if true: shift by shortest distance
    
    #if right of 1 is to left of left of other
    if r1[0] < r2[2] or r2[0] < r1[2]:
        return 0
    #if bottom of 1 is above top of other
    if r1[1] < r2[3] or r2[1] < r1[3]:
        return 0

    #shift vertically by 10% more than smallest vertical overlap
    v1 = r1[1] - r2[3]      
    v2 = r1[3] - r2[1]
    if shift_shortest:
        if abs(v1) < abs(v2):
            return -fac*v1
        else:
            return -fac*v2
    else:           #always shift vertically away from zero (origin of f-diag)
        #find y-loc of center of r1
        ctr_y = (r1[1] + r1[3])/2
        #if positive, shift positive, if negative shift negative
        #pick the dir that anti-aligns with ctr-y (diff sign)
        if ctr_y*v1 < 0:    #v1 is the way to shift
            return -fac*v1
        else:               #v2 is the way to shift
            return -fac*v2
        
    
#Place system point at origin
ogn = vec(0,0,0)
comment_loc = ogn - vec(0, sys_name_scale,0)

print('Place point representing system\'s center of mass and title it')
sys_pt = sphere(pos=ogn, radius=r_sys, color=sys_color)
    
#Place system title
f_h = 36        #font for labels
fscript_h = 16  #font for labels
s_script=0.5*f_h
sub_script=-0.25*f_h
mag_script=0.666*fscript_h
#kluge, create a sphere at the point where the system name goes for zoom purposes
sphere(pos=ogn+vec(0,sys_name_scale,0),color=vec(0,0,0),opacity=0)
sys_name = label(text=sys, pos=ogn+vec(0,sys_name_scale,0), align='center',
            height=f_h,box=False,line=False,color=sys_color,font='monospace')


#place forces: tail at ogn+0.5*r_sys in direction
Forces = []
Force_labels=[]
Force_label_spheres=[]
sfx = ['st','nd','rd','th']
for i, f in enumerate(Fvals):
    scene.waitfor('click')
    f_loc = ogn+r_sys*hat(f)  #needed for aesthetics and checking for overlap
    if i > 0:
        for g in Forces:
            if f_loc.equals(g.pos):
                f_loc = g.pos + g.axis
                #shift so it starts at end of next vector
                #rechecking means it will move farther out: don't break early
                #closest in always will be placed closest to system point
                #continuing will not skip closer vectors
    Forces.append(arrow(pos=f_loc, axis=f, shaftwidth=shaft_w,
        headwidth=head_wl, headlength=head_wl, 
        color=the_colors[i%len(the_colors)]))   #these colors never run...out
        
    if i >= 3:
        j = 3
    else:
        j = i
    print('Place '+str(i+1)+sfx[j]+' force vector')
    scene.waitfor('click')
    

    print('Place '+str(i+1)+sfx[j]+' force label')

    shift_labely = vec(0,-i*4*r_sys,0)   #in font-space, add to yoffset
    
    f_label_loc = vec(sys_name_scale,0.9*sys_name_scale,0) + shift_labely
    Force_label_spheres.append(sphere(pos=f_label_loc,radius=r_sys,
        color=vec(0,0,0),opacity=0))
    #preference to first vectors
    the_align = 'left'
    f_labelF = label(text='F', pos=f_label_loc, align=the_align,box=False,
                line=False, color=Forces[-1].color,height=f_h,font='monospace',
                opacity=0)

    f_labeltype = label(text=Fnames[i]['type'], pos=f_labelF.pos,
                    xoffset=s_script,yoffset=s_script, box=False, 
                    line=False,
                    align=the_align,color=Forces[-1].color,height=fscript_h,
                    font='monospace',opacity=0)

    f_labelobj = label(text=Fnames[i]['obj']+'->'+sys, pos=f_labelF.pos,
                    xoffset=s_script,yoffset=-s_script,box=False,
                    line=False,
                    align=the_align,color=Forces[-1].color,height=fscript_h,
                    font='monospace',opacity=0)

    print('Place magnitude, if known while solving (or add later once solved)')
    #print(f_h,f_labelF.height,f_labelobj.length,f_labelF.upper_right,f_labelF.lower_right,f_labelF.end, f_labelF.start)
    if show_magnitudes:
        fmag = mag(f)
        fmag_text = '='+str(fmag)+'N'
        if fmag > 20:    #only if its not too small a value
            if round(fmag) != fmag:
                fmag_text = '~'+str(round(fmag))+'N'
    
        f_labelmag = label(text=fmag_text,box=False,line=False,pos=f_labelobj.pos,
                    xoffset=len(f_labelobj.text)*mag_script,
                    yoffset=sub_script,
                    font='monospace',height=f_h,opacity=0,
                    align=the_align,color=Forces[-1].color)

    print('F^'+Fnames[i]['type']+'_'+Fnames[i]['obj']+'-->'+sys+'=',f,'N')
    if show_magnitudes:
        Force_labels.append([f_labelF,f_labeltype,f_labelobj,f_labelmag])
    else:
        Force_labels.append([f_labelF,f_labeltype,f_labelobj])

    
demo=label(text='click for tip-to-tail,tail-to-tip demo', 
        pos=ogn+vec(0,-sys_name_scale,0), align='center',
        height=f_h,box=False,line=False,color=sys_color,font='monospace')
scene.waitfor('click')
sys_name.visible=False
demo.visible=False
#for fl in Force_labels:
#    for l in fl:
#        l.visible=False

#Animate: tip-to-tail, tail-to-tip method (t6)
print('Duplicate forces for tip-to-tail, tail-to-tip method')
#simple estimate of max possible farthest right forces, scale
fmaxmag = 0
for f in Fvals:
    fmaxmag += mag(f)
#shift right by 3*fmaxmag

#duplicate forces (clone in same place then move original), no labels
t6Forces = []
for f in Forces:
    t6Forces.append(f.clone())
the_vel = vec(-1,0,0)

#Parameterize for animation:
dt = fmaxmag/shift_fac

scene.camera.follow(t6Forces[0])
scene.autoscale = False
while Forces[0].pos.x > ogn.x-3*fmaxmag:
    rate(shift_rate)
    #move them to the right for tip-to-tail
    for f in Forces:
        f.pos = f.pos + the_vel*dt
        
    sys_pt.pos = sys_pt.pos + the_vel*dt
        
    #exiting by position, no need to adjust time

scene.waitfor('click')

print('Move each force tip-to-tail onto previous. Shift opposite arrows slightly for clarity.')
#perform the method
t6shift = vec(0,0,0)    #store for later, when drawing Fnet, not perfect solution
#because if shifted multiple times, it could look funky but that's okay
for i,f in enumerate(t6Forces):
    if i > 0:
        #calculate how much time to take
        t = 0
        dt = 1
        f_last = t6Forces[i-1]
        f_last_tip = f_last.pos + f_last.axis
        f.v = (f_last_tip - f.pos)/shift_fac       #order matters, move position to last tip
        t_done = shift_fac
        #shift to side slightly if opposite directions
        ang_vec = abs(degrees(diff_angle(f.axis,f_last.axis)))
        if ang_vec > 160:
            #shift by average of shaft and headwidth
            #crossed with -z for counterclockwise shift
            shift_d = (9*shaft_w + head_wl)/10
            the_shift = shift_d*get_shift_hat(f_last.axis,False)
            if dot(the_shift,t6Forces[-1].pos+t6Forces[-1].axis) < 0:
                the_shift = -the_shift  #shift most aligned with radially out
            t6shift = t6shift + the_shift
            shift_v = the_shift/shift_fac
            t_done = shift_fac
            while t < t_done:
                rate(shift_rate)
                f.pos = f.pos + shift_v*dt
                t += dt
            t = 0   #reset for next
        while t < t_done:
            rate(shift_rate)
            f.pos = f.pos + f.v*dt
            t += dt

scene.waitfor('click')

print('Net force computed tail-to-tip: extends from tail of first to tip of last force')
print('This is shown as a dot if Fnet = 0, within computer error.')
#get the Fnet
Fnetval = Fvals[0]
for f in Fvals[1:]:
    Fnetval = Fnetval + f
Fnetmag = mag(Fnetval)
#allow for calculation error
feps = fmin*1e-5
fnet_align = 'left'
if Fnetmag < feps:
    Fnet = sphere(pos=t6Forces[0].pos+0.5*t6shift,radius=0.5*head_wl,color=net_color)
    Fnetmag = 0
    Fnet.axis = head_wl*hat(cross(t6Forces[-1].axis,vec(0,0,1)))
else:
    #Extend Fnet from t6Forces[0].pos to tip of last
    dt = 1
    t = dt
    t_done = shift_fac
    tip_last = t6Forces[-1].pos+t6Forces[-1].axis
    grow_v = (tip_last - t6Forces[0].pos - t6shift)/shift_fac

    ang_vec = abs(degrees(diff_angle(grow_v,t6Forces[-1].axis)))
    if ang_vec < 20:
        #shift by average of shaft and headwidth
        #crossed with -z for counterclockwise shift
        shift_d = (9*shaft_w+head_wl)/10
        the_shift = shift_d*get_shift_hat(t6Forces[-1].axis,False)
        if dot(the_shift,t6Forces[-1].pos+t6Forces[-1].axis) < 0:
            the_shift = -the_shift
        t6shift = t6shift + the_shift

    Fnet = arrow(pos=t6Forces[0].pos+t6shift, axis = grow_v*dt, shaftwidth=shaft_w, 
        headwidth=head_wl, headlength=head_wl,color=net_color)
    while t < t_done:
        rate(shift_rate)
        Fnet.axis = Fnet.axis + grow_v*dt
        
        t += dt

Fnet_label_shift = vec(0,0.9*sys_name_scale,0)

scene.autoscale=False
print('Label Fnet')
#label the Fnet

fnet_label_loc = Fnet.pos + Fnet_label_shift

fmag_text = 'Fnet'

if show_magnitudes:
    if Fnetmag > 20:    #only if its not too small a value
        if round(Fnetmag) != Fnetmag:
            fmag_text = fmag_text+'~'+str(round(Fnetmag))+'N'
        else:
            fmag_text = fmag_text+'='+str(Fnetmag)+'N'
    else:
        fmag_text = fmag_text+'='+str(Fnetmag)+'N'
        
print('Fnet =',Fnetval,'N')

the_align = 'left'
fnet_labelF = label(text=fmag_text, pos=fnet_label_loc, align='center',
            height=f_h,font='monospace',box=False,line=False,color=net_color,
            opacity=0)

