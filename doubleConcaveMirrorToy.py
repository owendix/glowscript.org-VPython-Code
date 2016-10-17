GlowScript 2.1 VPython
#demonstrates approximately what happens for the double concave mirror toy,
#where an image is placed inside the clam-shell of two concave mirrors, 
#the top with a hole in it, and the image of the object is ultimately
#projected to roughly where the hole is
#not using the thin optic approximation yields slightly better results, 
#though these still seem slightly high. Could 

thin_optic_approx = False

the_rate = 20  #for approx = false

silver = vec(0.75,0.75,0.75)
bottom_color = color.yellow
top_color = color.green
photon_opac = 0.5       #if thin_optic_approx = False
optics = []
h_o = 2
h_o2 = 0.666*h_o
r_o = 0.1*h_o
os_fac = 1.5
ax_fac = 2
c = 14
#make other arrow part half the height
obj_color1 = color.magenta
obj_color2 = color.cyan
img1_color = color.top_color
img2_color = color.bottom_color
img_opac = 0.8      #for False


canvw = 800
canvh = canvw
if thin_optic_approx:
    the_range = 1.5*c
else:
    the_range = 0.85*c
scn = canvas(width=canvw, height=canvh, forward=vec(0,0,-1), up=vec(0,1,0),
        autozoom=True,center=vec(0,0.1*c,0),range=the_range)



def make_arc(cs, h, sym_axis, sym_pt):
    global optics
    #make an arc from either 1 or 2 chords and the sagitta, or measurable height, 
    #if len(cs) == 1: chord and sagitta are used to easily make the arc
    #if len(cs) == 2: cs[0] > cs[1], the second chord says where to stop drawing
    #the arc, symmetrically about the sym_axis = vector from center of curvature
    #to symmetric center of arc, ASSUMES +Z IS OUT OF PAGE
    #sym_pt is the center of the arc, on the symmetry line
    #both follow from pythagorean theorem
    #Does not check for valid inputs
    sym_axis = hat(sym_axis)
    npts = 26       #must be even for symmetry, don't go too high, 26 is good
    #if npts%2 == 1:
    #    npts += 1       #even number of points
    if len(cs) == 1:
        c = cs[0]   #chord, diameter of part of sphere/circle
        c2 = c/2
        s = h       #sagitta is h, apothem is r - s
        #calculate radius of curvature
        r = 0.5*(c2**2/s + s)    #radius of chord, focus always r/2
        
        rsymvec = r*sym_axis
        rcvec = cross(vec(0,0,1),rsymvec)   #must be perpendicular, z+ out of page
        rpt = sym_pt - rsymvec
        
        #create 
        crvlist = []
        tmp=[0 for x in range(int(npts/2))]
        
        ang_max = asin(c/(2*r))
        dang = 2*ang_max/npts  

        #compute symmetric along symmetry line for increased accuracy
        for da in [dang, -dang]:    #must start out positive
            ang = da
            i = 0
            while abs(ang) <= ang_max:
                tmp[i]=rpt + cos(ang)*rsymvec + sin(ang)*rcvec
                ang = (i+2)*da      #less error than adding
                i += 1
            
            if da*dang > 0:         #same sign
                ang1 = ang
                tmp.reverse()       #doing this inside append works in glowscript, not python
                crvlist.extend(tmp)
                crvlist.append(sym_pt)
            else:
                crvlist.extend(tmp)
                
        optic_pieces = [curve(pos=crvlist,radius=c/200,color=bottom_color)]
        #append useful defining vectors, and the positive range of angles
        optic_pieces[0].params = [rpt, rsymvec, rcvec, 0, ang_max]
        #these ang limits are inclusive
        #sym_pt = rpt + rsymvec
    elif len(cs) == 2:
        c = cs[0]
        cstop = cs[1]
        #h = h
        r = sqrt((2*h**2 + c**2 + cstop**2)/8 + (c**4 + cstop**4 - 2*(c*cstop)**2)/(64*h**2))
        t = r - sqrt(r**2 - cstop**2/4)
        sym_pt = sym_pt + t*sym_axis
        
        rsymvec = r*sym_axis
        rcvec = cross(vec(0,0,1),rsymvec)   #must be perpendicular, z+ out of page
        rpt = sym_pt - rsymvec
        
        ang_min = asin(cstop/(2*r))
        ang_max = asin(c/(2*r))   #will be positive
        
        
        #do scale for same density of points
        npts = round((ang_max - ang_min)*npts/ang_max)
        if npts%2 == 1:
            npts += 1
        
        dang = 2*(ang_max-ang_min)/npts
            
        #create 
        crvlist = [0 for x in range(npts/2)]
        
        ang = ang_min
        #compute symmetric along symmetry line for increased accuracy
        optic_pieces = []
        new_ang_max = ang_max
        for da in [dang, -dang]:    #must start out positive
            i = 0
            while abs(ang) < (ang_max + dang):
                if ang > new_ang_max:
                    new_ang_max = ang
                crvlist[i] = rpt + cos(ang)*rsymvec + sin(ang)*rcvec
                ang += da      #less error than adding
                i += 1
            
            #if da*dang > 0:         #same sign
            #    crvlist.reverse()       #doing this inside append works in glowscript, not python
                
            ang = -ang_min
            optic_pieces.append(curve(pos=crvlist,radius=c/200,color=top_color))
            optic_pieces[-1].params=[rpt, rsymvec, rcvec, ang_min, new_ang_max]
            #ang limits are inclusive
            
    optics.append(optic_pieces)    #second optic in two pieces, produces an array of pieces
    
    return r

#make optics: list of lists of optic pieces [[bottom],[top_left],[top_right]]
#c = 14
#bottom mirror
s = 2.7
b_pt = vec(0,-4.76,0)
sym_axis = vec(0,-1,0)
rb = make_arc([c], s, sym_axis, b_pt)

#top mirror
t_pt = vec(0,0,0)
cstop=4.08
h = 1.99
sym_axis = -sym_axis        #symmetry axis points up now
rt = make_arc([c,cstop], h, sym_axis, t_pt)
#shift to symmetry point    #THIS EQUATION MAY BE WRONG, SHOULD SHIFT UP
t_pt = optics[-1][0].params[0] + optics[-1][0].params[1]

#draw axes:
if rb > rt:
    ax_lim = ax_fac*rb
else:
    ax_lim = ax_fac*rt
ax_r = ax_lim/500

axis = cylinder(pos=scn.center + vec(0,-ax_lim,0),axis=vec(0,2*ax_lim,0),radius=ax_r,
        opacity=0.5,color=color.white)
hlines = []
for pt in [b_pt, t_pt]:
    hlines.append(cylinder(pos=pt + vec(-ax_lim,0,0),axis=vec(2*ax_lim,0,0),
        radius=ax_r,opacity=0.25,color=color.white))
        
#put focii down
fs = [sphere(pos=b_pt+vec(0,rb/2,0),radius=3*ax_r,color=bottom_color),
        sphere(pos=t_pt+vec(0,-rt/2,0),radius=3*ax_r,color=top_color)]
fs_label = [label(text='bottom mirror focus',color=bottom_color),
            label(text='top mirror focus',color=top_color)]
for i, l in enumerate(fs_label):
    l.pos = fs[i].pos + vec(-c/2,0,0)
    l.align='right'
    l.box=False
    l.line=False
    l.height=18
    l.opacity=0
    
#put object down
obj_shift = os_fac*t_pt.y

#sym_axis is up
obj_pos = b_pt + obj_shift*sym_axis + vec(-h_o/2,0,0)
object=[cylinder(pos=obj_pos, axis=vec(h_o,0,0),radius=r_o,color=obj_color1),
        cylinder(pos=obj_pos, axis=vec(0,0,h_o2),radius=r_o,color=obj_color2)]


if thin_optic_approx:
    #make first image, using top mirror, distance from t_pt.y along y
    approx_disclaimer = label(text='Uses the thin optic approximation',
        height=24, pos=scn.center-vec(0,0.7*scn.range,0), align='center', 
        box=False, line=False,opacity=0, color=color.red)
    o1 = abs(obj_pos.y - t_pt.y)
    f1 = rt/2
    
    i1 = o1*f1/(o1 - f1)
    io = -i1/o1
    h_i = io*object[0].axis.x
    h_itail = io*object[0].pos.x
    h_i2 = h_i/2
    
    r_i1 = 0.1*h_i
    i1pos = t_pt + vec(h_itail,-i1,0)
    image1 = [cylinder(pos=i1pos,axis=vec(h_i,0,0),radius=r_i1,color=obj_color1,
                opacity=0.333),
                cylinder(pos=i1pos,axis=vec(0,0,h_i2),radius=r_i1,color=obj_color2,
                opacity=0.333)]
    
    #make rays for first image
    o_tip = object[0].pos + object[0].axis
    impact = vec(o_tip.x, hlines[1].pos.y,0)
    rays1a=[cylinder(pos=o_tip,axis=impact-o_tip,radius=ax_r,
            color=top_color)]
    rays1a.append(cylinder(pos=impact,axis=2*(fs[1].pos - impact),radius=ax_r,
            color=top_color))
    i_tip = image1[0].pos + image1[0].axis
    rays1a.append(cylinder(pos=impact,axis=1.2*(i_tip - impact),radius=ax_r,
            color=top_color,opacity=0.5))
    impact = vec(axis.pos.x,hlines[1].pos.y,0)
    rays1b = [cylinder(pos=o_tip, axis=impact - o_tip, radius = ax_r,
                color=top_color)]
    rays1b.append(cylinder(pos=impact, axis=2*(vec(-o_tip.x,o_tip.y,0)-impact),radius=ax_r,
            color=top_color))
    rays1b.append(cylinder(pos=impact,axis=1.2*(i_tip - impact),radius=ax_r,
            color=top_color,opacity=0.5))
    
    #make second image, using bottom mirror
    o2 = abs(image1[0].pos.y - b_pt.y)
    f2 = rb/2
    
    i2 = o2*f2/(o2 - f2)
    io = -i2/o2
    h_i = h_i = io*image1[0].axis.x
    h_itail = io*image1[0].pos.x
    h_i2 = h_i/2
    r_i2 = 0.1*h_i
    i2pos = b_pt + vec(h_itail,i2,0)
    image2 = [cylinder(pos=i2pos,axis=vec(h_i,0,0),radius=r_i2,color=obj_color1,
                opacity=0.666),
                cylinder(pos=i2pos,axis=vec(0,0,h_i2),radius=r_i2,color=obj_color2,
                opacity=0.666)]
    
    #make rays for second image
    o_tip = image1[0].pos + image1[0].axis
    impact = vec(o_tip.x, hlines[0].pos.y,0)
    rays2a=[cylinder(pos=o_tip,axis=impact-o_tip,radius=ax_r,
            color=bottom_color)]
    rays2a.append(cylinder(pos=impact,axis=2*(fs[0].pos - impact),radius=ax_r,
            color=bottom_color))
    impact = vec(axis.pos.x,hlines[0].pos.y,0)
    rays2b = [cylinder(pos=o_tip, axis=impact - o_tip, radius = ax_r,
                color=bottom_color)]
    rays2b.append(cylinder(pos=impact, axis=2*(vec(-o_tip.x,o_tip.y,0)-impact),radius=ax_r,
            color=bottom_color))
    
    show_rays = [True, True]
    
    def Show_Rays(bray):
        global show_rays, rays1a, rays1b, rays2a, rays2b
        
        if '1' == bray.text[-1]:
            if 'Hide' in bray.text:
                show_rays[0] = False
                bray.text = 'Show Rays for Image1'
                for r in rays1a:
                    r.visible = False
                for r in rays1b:
                    r.visible = False
                    
            else: #show
                show_rays[0] = True
                bray.text = 'Hide Rays for Image1'
                for r in rays1a:
                    r.visible = True
                for r in rays1b:
                    r.visible = True
        else:   #2nd
            if 'Hide' in bray.text:
                show_rays[1] = False
                bray.text = 'Show Rays for Image2'
                for r in rays2a:
                    r.visible = False
                for r in rays2b:
                    r.visible = False
                    
            else: #show
                show_rays[1] = True
                bray.text = 'Hide Rays for Image2'
                for r in rays2a:
                    r.visible = True
                for r in rays2b:
                    r.visible = True
    
    scn.title = 'Toggle Principle Rays On/Off\n\n'
    #button ray1
    bray1 = button(text='Hide Rays for Image1', pos=scn.title_anchor, bind=Show_Rays)
    #button ray1
    bray2 = button(text='Hide Rays for Image2', pos=scn.title_anchor, bind=Show_Rays)
    scn.append_to_title('\n')
    
            
else:   
    #not the thin optic approximation
    dt = 2**-10      #power of 2 good for computer accuracy
    mirror_frac = 0.7           #=0.7, with npts=26 in mirror, don't change, its good
    n_photons = int(2*optics[1][0].npoints*mirror_frac)-1
    #start at the tip of the object
    o_tip = object[0].pos + object[0].axis
    photons = [0 for q in range(n_photons)]
    n = 0
    while n < n_photons:
        rate(the_rate)
        #pick a direction that's between the beginning and end of the 
        #top optic: optics[1][0] and optics[1][1]
        #ref. bottom optic: optics[0][0]
        rphoton = o_tip
        q = vec.random()
        #pick a random curve in curvelist
        side = round(abs(q.x))
        the_optic = optics[1][side]

        #which point to aim at first
        #print('side:',side,the_optic.npoints)
        i_mirror = int(0.8*the_optic.npoints*abs(q.y))   #aim at half nearest to sym_pt
        impact_pt = the_optic.point(i_mirror).pos   #where it will hit
        the_dir = impact_pt - o_tip
        if n > 0:
            same_dir = True
            try_again = False
            while same_dir:
                for p in photons[:n]:
                    if the_dir.equals(p[0].axis):
                        #shot at same direction, pick a different one
                        q=vec.random()
                        side = round(abs(q.x))
                        the_optic = optics[1][side]
                        i_mirror = int(0.5*the_optic.npoints*abs(q.y))
                        #print('same dir:',i_mirror)
                        impact_pt = the_optic.point(i_mirror).pos   #where it will hit
                        the_dir = impact_pt - o_tip
                        try_again = True
                        break
                if not try_again:
                    same_dir = False    #found a unique direction
                else:
                    try_again = False   #reset for new direction
                
                
        photons[n] = [cylinder(pos=o_tip, axis=the_dir, 
                    radius=ax_r, color=obj_color1,opacity=photon_opac)]
        v = hat(photons[n][0].axis)
        [rpt, rsym, rc, ang_min, ang_max] = the_optic.params
        #once first photon created as a cylinder, it has hit the top mirror

        #determines when/if it will bounce
        last_op = 1
        in_mirror = True
        n_bounce = 0
        while in_mirror:
            rate(the_rate)
            n_bounce += 1
            v = v - 2*proj(v,impact_pt - rpt)
            v = hat(v)  #in case of computer error
            rphoton = impact_pt + v*dt    #to avoid t = 0 as only solution
            #while loop to get there, once it gets there enter another loop-set
            #that breaks when the photon leaves the cavity and make it move a little
            #past that
            #see which optic this new beam will bounce off of
            #top optic hit last, try bottom first, if bottom last, top first
            step = 2*last_op - 1    # 1 or -1
            start = (last_op + 1)%2
            #perhaps try using the numerical method
            #imagine a line segment of light, intersecting little line segments
            #still draw a radius to that point, once impacted to get how it
            #bounces
            for i_op in range(start, last_op + step, step):
                [rpt, rsym, rc, ang_min, ang_max] = optics[i_op][0].params

                r2 = mag2(rsym)
                rph_pt = rphoton-rpt
                C = mag2(rph_pt) - r2    #always negative if inside mirror's full circle
                if C > 0:   #I shifted slightly to avoid C=0, same optic
                    continue    #probably unnecessary
                #optics[#][#].params = [rpt, rsymvec, rcvec, ang_min, ang_max] #angs inclusive
                #dot(rsymvec,rcvec)=0, mag(rsymvec) = mag(rcvec) = r
                #symmetry point in optic is rpt + rsymvec
                            
                Bhalf = dot(rph_pt,v)
                #A = 1 (mag2(v) = 1)
                
                t = -Bhalf + sqrt(Bhalf**2 - C) #only pos one will be non-neg
                
                #use sine term to find angle, it has the correct domain
                rtv = rph_pt + t*v
                c_ang = dot(rtv,rsym)/r2        #this is not always less than 1!!!
                if c_ang < 0:   #definitely aimed at other one
                    continue    #if false, still may be other one
                s_ang = dot(rtv,rc)/r2
                ang = asin(s_ang)
                
                if i_op == 1 and ang < ang_min and ang > -ang_min:
                    #it will bounce out
                    in_mirror = False
                    #project forward a bit
                    photons[n].append(cylinder(pos=rphoton,axis=sqrt(r2)*v,
                                    radius=ax_r, color=obj_color1,
                                    opacity=photon_opac))
                    break
                elif ang < ang_max and ang > -ang_max:  #it will bounce
                    impact_pt = rpt + c_ang*rsym + s_ang*rc
                    photons[n].append(cylinder(pos=rphoton,axis=(impact_pt-rphoton),
                            radius=ax_r, color=obj_color1,opacity=photon_opac))
                    last_op = i_op
                    break
                else:
                    #it won't bounce here
                    continue
                
        n += 1  #go to next photon
    
    #got all the photons, photons[0-n_photons-1][-1] is last beam exiting
    r_isct=[]
    for i, p in enumerate(photons[:-1]):
        #find where they intersect each other
        for q in photons[i+1:]:
            rp = p[-1].pos
            ap = p[-1].axis
            rq = q[-1].pos
            aq = q[-1].axis
            xqp_pos = rq.x - rp.x
            yqp_pos = rq.y - rp.y
            xypa = ap.x/ap.y
            t_isct = (xqp_pos - yqp_pos*xypa)/(aq.y*xypa - aq.x)
            r_isct.append(rq + t_isct*aq)
    r_avg = vec(0,0,0)
    for r in r_isct:
        r_avg = r_avg + r
    #not a numerically stable method but lets run with it
    r_avg = r_avg/len(r_isct)
    click_l = label(text='click to show image at average convergence point',
        height=24, pos=scn.center-vec(0,0.7*scn.range,0), align='center', 
        box=False, line=False,opacity=0, color=color.red)
    scn.waitfor('click')
    click_l.text='image is at average convergence point'
    sphere(pos=r_avg,radius=3*ax_r,color=color.red)   #just to mark avg convergence point
    #2*r_avg.x is the height
    i_tip = r_avg
    i_axis = vec(2*r_avg.x,0,0)
    i_pos = i_tip - i_axis
    h_i = i_axis.x      #if this is negative and other is positive, will flip
    h_ifac = h_i/object[0].axis.x
    h_i2 = h_ifac*object[1].axis.z
    r_i = abs(h_ifac)*object[0].radius
    the_image = [cylinder(pos=i_pos,axis=i_axis,radius=r_i,color=obj_color1,
            opacity=img_opac),
            cylinder(pos=i_pos,axis=vec(0,0,h_i2),radius=r_i,color=obj_color2,
            opacity=img_opac)]
