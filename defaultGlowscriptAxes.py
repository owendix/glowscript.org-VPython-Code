GlowScript 2.1 VPython
#Show the 3-axis for the default coordinate system glowscript makes
#Use different colors for x,y,z and label + and -
#Able to include fewer axes, and just positive or negative directions
#Can also add a ball at a point with different velocity, accelerations


#Adding a ball
ball_pos = vec(1,2,0)

put_ball = False
put_pos_arrow = False
put_vel_arrow = False
put_acc_arrow = False

if put_vel_arrow:
    ball_vel = vec(1.5,0,0)
if put_acc_arrow:
    ball_acc = vec(0,-1,0)
sw = 0.05
ball_rad = 0.1
ball_color=color.yellow
eps_fact = 1.15
#end adding a ball code


#Adjusting the axes
#uses this origin in seting axes, other features (like axis length) not controllable
the_origin=vec(0,0,0)  #use this
the_lim = 4
the_make_xyz = [True,True,True]     #make all 3 axes if all true
the_make_pn = [True,True]           #make positive & negative dirs if both true
#end adjusting axes


#controling view
canv_width = 1200
canv_height = 9*canv_width/16
#don't adjust scene to origin, may want to show difference
scene = display(width=canv_width,height=canv_height)


#picking values gives default overrideable when called farther down page
def putAxes(origin, lim, make_xyz, make_pn):
    """                
    Makes basic axes, lim units long with (int(lim)-1) tick marks
    """
    
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
                shaftwidth=sw,headwidth=hwl,headlength=hwl))
            for tf in tick_factors:
                #apply tick marks (thin square solids)
                the_ticks.append(box(pos=origin+tf*pn*dir,axis=pn*dir,length=tL,
                    height=tHW,width=tHW,color=color.white))
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

            
putAxes(the_origin, the_lim, the_make_xyz, the_make_pn)

if put_ball:
    ball=sphere(pos=ball_pos,radius=ball_rad, color=ball_color)

    if put_pos_arrow:
        ball_pos_vec=arrow(pos=the_origin,axis=ball.pos,shaftwidth=sw,
            headwidth=1.5*ball_rad,headlength=2*ball_rad,color=color.red)
    if put_vel_arrow:
        ball_vel_vec=arrow(pos=ball.pos,axis=ball_vel,shaftwidth=sw,
            color=color.blue)
        text(text='vel',pos=ball.pos+ball_vel_vec.axis*eps_fact,
                height=2*ball_rad, align='center',color=color.blue)
        
    if put_acc_arrow:
        ball_acc_vec = arrow(pos=ball.pos,axis=ball_acc,shaftwidth=sw,
            color=color.white)
        text(text='acc',pos=ball.pos+ball_acc_vec.axis*eps_fact,
            height=2*ball_rad, align='center',color=color.white)
