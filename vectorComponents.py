GlowScript 2.1 VPython
#Show the 3-axis for the default coordinate system glowscript makes
#Use different colors for x,y,z and label + and -
#Able to include fewer axes, and just positive or negative directions
#Can also add a ball at a point with different velocity, accelerations

vec_axis = vec(2,5,0)
the_lim = 7

comp_order = 'xyz'      #'xyz', 'xzy', 'yzx', 'yxz', 'zxy', 'zyx'

#Adjusting the axes
#uses this origin in seting axes, other features (like axis length) not controllable
the_make_xyz = [True,True,True]     #make all 3 axes if all true
the_make_pn = [True,True]           #make positive & negative dirs if both true
#end adjusting axes

shaft_w = 0.05*the_lim
head_wl = 2*shaft_w
vec_text_h = head_wl
grow_rate = 1000
grow_fac = 50

#controling view
ogn=vec(0,0,0)  #use this
canv_width = 800
canv_height = 9*canv_width/16
#don't adjust scene to origin, may want to show difference
scene = display(width=canv_width,height=canv_height,background=color.black)


#picking values gives default overrideable when called farther down page
def putAxes(origin, lim, make_xyz, make_pn):
    """                
    Makes basic axes, lim units long with (int(lim)-1) tick marks
    """
    opac = 0.2
    lim_min = 2
    lim_max = 10
    if lim < lim_min:
        lim = lim_min
    elif lim > lim_max:
        lim = lim_max
        
    dirs = []                    #=[lim*vec(1,0,0),lim*vec(0,1,0),lim*vec(0,0,1)]
    pns = []                     #=[1,-1]
    dirs_str = []                #=['x','y','z']
    pns_str = []                 #=['+','-']
    if make_xyz[0]:
        dirs.append(lim*vec(1,0,0))
        dirs_str.append('x')
    if make_xyz[1]:
        dirs.append(lim*vec(0,1,0))
        dirs_str.append('y')
    if make_xyz[2]:
        dirs.append(lim*vec(0,0,1))
        dirs_str.append('z')
    if make_pn[0]:
        pns.append(1)
        pns_str.append('+')
    if make_pn[1]:
        pns.append(-1)
        pns_str.append('-')        
    
    #print('dir, pn, dir_str, pn_str: ',dir, pn, dir_str, pn_str)
    
    #Don't really need to store these lists, but in case I want to use them
    the_axes=[]
    the_labels=[]
    the_ticks=[]
    tick_factors=[]
    for i in range(1,int(lim)):
        tick_factors.append(i/int(lim))     #=[1/3, 2/3]     factors
    
    sw = 0.0333*lim     #axes shaftwidth
    hwl = 2*sw          #axes headwidth, headlength same
    tL = 0.2*sw
    tHW = 2*sw
    ht_label = 0.1*lim
    eps_fact = 1.1  #put text out a little past end of axes
    

    for i_dir, dir in enumerate(dirs):
        for i_pn, pn in enumerate(pns):
            #put the axes arrows
            the_axes.append(arrow(pos=origin,axis=pn*dir,color=dir,
                shaftwidth=sw,headwidth=hwl,headlength=hwl, opacity=opac))
            for tf in tick_factors:
                #apply tick marks (thin square solids)
                the_ticks.append(box(pos=origin+tf*pn*dir,axis=pn*dir,length=tL,
                    height=tHW,width=tHW,color=color.white,opacity=opac))
                #No tick labels included
            #access last one as the_axes[-1]
            #labels the axes
            the_labels.append(text(pos=origin+eps_fact*pn*dir,
                text=pns_str[i_pn]+dirs_str[i_dir],align='center',height=ht_label,
                color=dir, billboard=True)) #billboard makes text face camera

    if len(dirs) == 1:
        #if just one dimension, include an origin tick, labeled as 0
        ogn=box(pos=origin,axis=dir,length=tL,height=tHW,width=tHW,
            color=color.white)
        #if y or z, shift right, if x shift down
        if dot(dirs[0],vec(1,0,0)) != 0:    #points in x-dir
            shift=vec(0,-2,0)
        else:
            shift=vec(1,0,0)
        text(pos=origin+shift*tHW*eps_fact,text='0',height=ht_label,
            align='center',billboard=True,color=color.white)
    
    return      #not needed but emphasizes function is done and returns nothing
            
putAxes(ogn, the_lim, the_make_xyz, the_make_pn)

the_vec = arrow(pos=ogn,axis=vec_axis,shaftwidth=shaft_w, headwidth=head_wl,
                headlength=head_wl,color=color.yellow)

vec_text = '('+vec_axis.x+','+vec_axis.y+','+vec_axis.z+')'
print('Vector: '+vec_text)
comps = [the_vec]

vec_loc = ogn
eps = 1e-2
dir_dict = {'x': vec(1,0,0), 'y': vec(0,1,0), 'z': vec(0,0,1)}
dir_text = list(comp_order)     #e.g. ['x', 'y', 'z']
dirs = [dir_dict[d] for d in dir_text]

for i, dir in enumerate(dirs):
    comp_axis = proj(vec_axis,dir)
    if mag(comp_axis) > eps:
        scene.waitfor('click')
        dt = 1
        t = dt
        comp_v = comp_axis/grow_fac
        t_done = grow_fac
        comp_text = '('+comp_axis.x+','+comp_axis.y+','+comp_axis.z+')'
        print(dir_text[i]+'-component vector: '+comp_text)
        comps.append(arrow(pos=vec_loc, axis=comp_v*dt, shaftwidth=shaft_w, 
                    headwidth=head_wl, headlength=head_wl, color=dir))
        while t < t_done:
            rate(grow_rate)
            comps[-1].axis = comps[-1].axis + comp_v*dt
            t += dt            

        vec_loc = vec_loc + comp_axis
    
                
                
