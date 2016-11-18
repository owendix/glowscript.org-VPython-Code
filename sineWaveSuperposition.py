GlowScript 2.2 VPython
#shows the superposition of two sine waves, one with adjustable amplitude, 
#period, and phase shift

the_rate0 = 100 #
the_rate = 1.0*the_rate0

Nwaves = 4      #can include superposition of many waves
#min is 2 waves, more than 5 and it slows down a fair amount
if Nwaves < 2:
    Nwaves = 2
elif Nwaves > 5:
    Nwaves = 5
wave_radfac = 1/800
sumwave_radfac = 3
sumwave_color = color.white
slide_fac = 4
slide_step = 2e-4
font_h = 18
big_sldr_lenfac = 3.0075

y_max = 1
x_max = 4


cw = 1100
ch = 0.65*cw
scn = display(width=cw,height=ch,center=vec(0.4*x_max,0,0),range=0.6*x_max,
            userzoom=False,userspin=False)
scn.autoscale=False

#create axes
axes_color = color.gray(0.75)
axes_sw = 1.5*wave_radfac*x_max
axes_hl = 10*axes_sw
axes_hw = 0.5*axes_hl
xax_pad = 1.05
yax_pad = 1.5
ogn = vec(0,0,0)
ax_xp = arrow(pos=ogn,axis=xax_pad*vec(x_max,0,0),shaftwidth=axes_sw,
            headlength=axes_hl,headwidth=axes_hw,color=axes_color)
ax_xlabel = label(text='x(m)',pos=ogn+ax_xp.axis,yoffset=-font_h,align='center',
            box=False,line=False,opacity=0,height=font_h,color=axes_color)
ax_yp = arrow(pos=ogn,axis=yax_pad*vec(0,y_max,0),shaftwidth=axes_sw,
            headlength=axes_hl,headwidth=axes_hw,color=axes_color)
ax_ylabel = label(text='y',pos=ogn+ax_yp.axis,yoffset=0.5*font_h,align='center',
            box=False,line=False,opacity=0,height=font_h,color=axes_color)
ax_yn = arrow(pos=ogn,axis=yax_pad*vec(0,-y_max,0),shaftwidth=axes_sw,
            headlength=axes_hl,headwidth=axes_hw,color=axes_color)


#initialize waves
t = 0   #global time variable
t_color = color.yellow
spacebar_txt = 's\nhit spacebar\nto (un)pause'
t_lbl = label(text='t='+'{:.2f}'.format(t)+spacebar_txt,pos=ogn+ax_xp.axis+ax_yp.axis,
            align='left',box=False,line=False,opacity=0,height=font_h,
            color=t_color)
dt = 5e-3
dx = 8e-3
v = 0.5   #all waves traveling same speed, allow opposite directions: L = vT
waves = [0,0]
A0s = [1,1]
T0s = [2,2]
shift0s = [0,pi/2]
dir0s = [1,1]
wave_colors0 = [color.cyan, color.yellow, color.magenta, color.green, color.red, color.orange, color.blue]
wave_colors = wave_colors0[0:2].copy()
for n in range(2,Nwaves):
    waves.append(0)
    A0s.append(0)   #make the amplitudes 0 after first 2 waves
    T0s.append(2)
    shift0s.append(2*pi*abs(vec.random().x))
    dir0s.append(1)
    wave_colors.append(wave_colors0[n%len(wave_colors0)])

#A0s = [1 for x in range(Nwaves)]
#T0s = [1 for x in range(Nwaves)]
#shift0s = [0,pi/2,0.6*pi,1.5*pi,pi/4,3*pi/4,5*pi/4,7*pi/4]    #
    
c_pos = []
pts_base = []

As = A0s.copy()
Ts = T0s.copy()
shifts = shift0s.copy()
dirs = dir0s.copy()



def change_wave(i,A,T,shift,dir):
    #assumes c_pos, pts_base, wave already exist
    global c_pos, waves
    x = 0
    n = 0
    while x < x_max:
        waves[i].modify(n,y=A*sin((2*pi/Ts[i])*(x/v - dirs[i]*t) - shifts[i]))
        
        n += 1
        x = n*dx    #more stable way to add

for i,A in enumerate(As):
    x = 0
    n = 0
    if i == 0:
        while x < x_max:
            c_pos.append(vec(x,A*sin((2*pi/Ts[i])*(x/v - dirs[i]*t) - shifts[i]),0))
            pts_base.append(vec(x,0,0))
            n += 1
            x = n*dx    #more stable way to add
    else:
        x = 0
        n = 0
        while x < x_max:
            c_pos[n] = vec(x,A*sin((2*pi/Ts[i])*(x/v - dirs[i]*t) - shifts[i]), 0)
            
            n += 1
            x = n*dx    #more stable way to add
            
    waves[i] = curve(c_pos)
    waves[i].radius = wave_radfac*x_max
    waves[i].color = wave_colors[i]
    if i >= 2:
        waves[i].visible=False
    
sumwave = 0

def init_sumwave():
    global pts_base, sumwave
    
    for i in range(len(pts_base)):
        pts_base[i].y = 0
        for w in waves:
            pts_base[i].y += w.point(i).pos.y
    
    sumwave = curve(pts_base)
    sumwave.radius=sumwave_radfac*wave_radfac*x_max
    sumwave.color = sumwave_color
            
init_sumwave()

def change_sumwave():
    global pts_base, sumwave
    
    for i in range(len(pts_base)):
        pts_base[i].y = 0
        for w in waves:
            pts_base[i].y += w.point(i).pos.y
            
        sumwave.modify(i,y=pts_base[i].y)

#add labels for waves
wave0_txt = 'wave1\nA='+str(As[0])+'\nL='+str(v*Ts[0])+'m\nT='+str(Ts[0])+'s\ns='+str(shifts[0])
wave0_txt += '\ndir='+str(dirs[0])
wave0_lbl = label(text=wave0_txt, pos=vec(0,y_max,0),box=False,line=False,
                opacity=0,height=font_h,
                xoffset=-5,align='right',font='monospace',color=wave_colors[0])


eqn_lbl = label(text='y=Asin(2pi x/L - dir(2pi t/T) - s)',
                pos=vec(0,y_max/5,0),box=False,line=False,opacity=0,height=(7/9)*font_h,
                xoffset=-5,align='right',font='monospace',color=color.white)
if Nwaves <= 2:
    sumwave_txt = '\nL = vT\nv='+'{:.1f}'.format(v)+'m/s\nytot = y1 + y2\n'
else:
    sumwave_txt =  '\nL = vT\nv='+'{:.1f}'.format(v)+'m/s\nytot = y1 + y2 + ...\n'
sumwave_lbl = label(text=sumwave_txt,pos=vec(0,-y_max/5,0),box=False,line=False,
                height=font_h,opacity=0,
                xoffset=-5,align='right',font='monospace',color=sumwave_color)

wave1_txt = 'wave2\nA='+'{:.2f}'.format(As[1])+'\nL='+'{:.2f}'.format(v*Ts[1])
wave1_txt += 'm\nT='+'{:.2f}'.format(Ts[1])+'s\ns='+'{:.2f}'.format(shifts[1])
wave1_txt += '\ndir='+str(dirs[1])
wave1_lbl = label(text=wave1_txt, pos=vec(0,-y_max,0),box=False,line=False,
                opacity=0,height=font_h,
                xoffset=-5,align='right',font='monospace',color=wave_colors[1])

#buttons for wave visibility
def toggle_wave_visibility(the_radio):    
    global waves
    
    i = int(the_radio.text[-1]) - 1

    waves[i].visible = not waves[i].visible
    if waves[i].visible:
        the_radio.text = 'Hide wave'+str(i+1)
    else:
        the_radio.text = 'Show wave'+str(i+1)
    
wave_vis = []

wave_vis.append(button(text='Hide wave1', pos=scn.title_anchor, 
            bind=toggle_wave_visibility))

wave_vis.append(button(text='Hide wave2', pos=scn.title_anchor, 
            bind=toggle_wave_visibility))

for i in range(2,Nwaves):
    wave_vis.append(button(text='Show wave'+str(i+1),pos=scn.title_anchor, 
            bind=toggle_wave_visibility))

            

def set_A(s):
    global As, wave1_lbl
    i = s.wavenum
    As[i] = s.value
    if i == 1:
        wave1_txt = 'wave2\nA='+'{:.2f}'.format(As[1])+'\nL='+'{:.2f}'.format(v*Ts[1])
        wave1_txt += 'm\nT='+'{:.2f}'.format(Ts[1])+'s\ns='+'{:.2f}'.format(shifts[1])
        wave1_txt += '\ndir='+str(dirs[1])
        wave1_lbl.text=wave1_txt
    
    change_wave(i,As[i],Ts[i],shifts[i],dirs[i])
    change_sumwave()
    
def set_T(s):
    global Ts, wave1_lbl
    i = s.wavenum
    Ts[i] = s.value
    if i == 1:
        wave1_txt = 'wave2\nA='+'{:.2f}'.format(As[1])+'\nL='+'{:.2f}'.format(v*Ts[1])
        wave1_txt += 'm\nT='+'{:.2f}'.format(Ts[1])+'s\ns='+'{:.2f}'.format(shifts[1])
        wave1_txt += '\ndir='+str(dirs[1])
        wave1_lbl.text=wave1_txt
        
    change_wave(i,As[i],Ts[i],shifts[i],dirs[i])
    change_sumwave()
    
def set_shift(s):
    global shifts, wave1_lbl
    i = s.wavenum
    shifts[i] = s.value
    if i == 1:
        wave1_txt = 'wave2\nA='+'{:.2f}'.format(As[1])+'\nL='+'{:.2f}'.format(v*Ts[1])
        wave1_txt += 'm\nT='+'{:.2f}'.format(Ts[1])+'s\ns='+'{:.2f}'.format(shifts[1])
        wave1_txt += '\ndir='+str(dirs[1])
        wave1_lbl.text = wave1_txt
        
    change_wave(i,As[i],Ts[i],shifts[i],dirs[i])
    change_sumwave()

def set_dir(r):
    global dirs, wave1_lbl
    
    i = r.wavenum
    #r.checked = not r.checked  #it switches this automatically when clicked
    #checked means wave travels in the negative direction
    if r.checked:
        dirs[i] = -1
    else:
        dirs[i] = 1
        
    if i == 1:
        wave1_txt = 'wave2\nA='+'{:.2f}'.format(As[1])+'\nL='+'{:.2f}'.format(v*Ts[1])
        wave1_txt += 'm\nT='+'{:.2f}'.format(Ts[1])+'s\ns='+'{:.2f}'.format(shifts[1])
        wave1_txt += '\ndir='+str(dirs[1])
        wave1_lbl.text = wave1_txt
            
    change_wave(i,As[i],Ts[i],shifts[i],dirs[i])
    change_sumwave()

#add slider to change animation rate:
def set_rate(s):
    global the_rate
    the_rate = s.value

scn.caption='Vary the animation rate:\n'

sldr_len = 0.275*cw
rate_fac = 20
slider(min=the_rate0/rate_fac,value=the_rate,max=the_rate0*rate_fac,
        length=big_sldr_lenfac*sldr_len,bind=set_rate)

scn.append_to_caption('\n')

#add sliders for second wave and beyond, first wave unshifted
A_sldrs = []
T_sldrs = []
shift_sldrs = []
dir_radios = []
#sldr_len = 0.275*cw
for i in range(1,Nwaves):
    scn.append_to_caption('Vary wave'+str(i+1)+' Amplitude: \t\t\t\t\t Vary wave')
    scn.append_to_caption(str(i+1)+' Period: \t\t\t\t\t Shift wave'+str(i+1)+':\n')
    
    #use A0s[0] for baseline min and max: A0s[2 or higher] start at 0
    A_sldrs.append(slider(min=0,value=A0s[i],max=A0s[0]*slide_fac,
                    pos=scn.caption_anchor,step=slide_step,length=sldr_len,bind=set_A))
    A_sldrs[i-1].wavenum = i
    #scn.append_to_caption(' ') #doesn't do anything
    
    #use T0s[0] for baseline min and max
    T_sldrs.append(slider(min=T0s[0]/(10*slide_fac),value=T0s[i],max=T0s[0]*slide_fac,
                    pos=scn.caption_anchor,step=slide_step,length=sldr_len,bind=set_T))
    T_sldrs[i-1].wavenum = i

    #scn.append_to_caption(' ') #doesn't do anything
    
    shift_sldrs.append(slider(min=0,value=shift0s[i],max=2*pi,
                    pos=scn.caption_anchor,length=sldr_len,bind=set_shift))
    shift_sldrs[i-1].wavenum = i
    
    scn.append_to_caption('\t')
    
    dir_radios.append(radio(checked=(dirs[i]==-1),text='wave'+str(i+1)+' neg dir?',bind=set_dir))
    dir_radios[i-1].wavenum = i

    scn.append_to_caption('\n')
    
#spacebar pauses all motion: set is_paused
is_paused = True

def pause_waves(evt):
    global is_paused
    s = evt.key
    if s == ' ':
        is_paused = not is_paused
        
#uses keyup so it can't be held down
scn.bind('keyup',pause_waves)


#animate waves
spacebar_txt = 's\nhit spacebar\nto (un)pause'
while True:
    rate(the_rate)
    if not is_paused:
        for i in range(len(waves)):
            change_wave(i,As[i],Ts[i],shifts[i],dirs[i])
        change_sumwave()
        t_lbl.text='t='+'{:.2f}'.format(t)+spacebar_txt
        t += dt

