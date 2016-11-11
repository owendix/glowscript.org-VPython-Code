GlowScript 2.1 VPython
#shows the superposition of two sine waves, one with adjustable amplitude, 
#period, and phase shift

wave_colors = [color.cyan, color.yellow]#, color.magenta]
wave_radfac = 1/800
sumwave_radfac = 3
sumwave_color = color.red
slide_fac = 4
slide_step = 2e-4
font_h = 18

y_max = 1
x_max = 4

Nwaves = 2      #only works for 2 right now

cw = 900
ch = 0.8*cw
scn = display(width=cw,height=ch,center=vec(0.5*x_max,0,0),range=0.65*x_max,
            userzoom=False,userspin=False)
scn.autoscale=True

#create axes
axes_color = color.white
axes_sw = 0.75*wave_radfac*x_max
axes_hl = 10*axes_sw
axes_hw = 0.5*axes_hl
axes_pad = 1.2
ogn = vec(0,0,0)
ax_xp = arrow(pos=ogn,axis=axes_pad*vec(x_max,0,0),shaftwidth=axes_sw,
            headlength=axes_hl,headwidth=axes_hw,color=axes_color)
ax_yp = arrow(pos=ogn,axis=axes_pad*vec(0,y_max,0),shaftwidth=axes_sw,
            headlength=axes_hl,headwidth=axes_hw,color=axes_color)
ax_yn = arrow(pos=ogn,axis=axes_pad*vec(0,-y_max,0),shaftwidth=axes_sw,
            headlength=axes_hl,headwidth=axes_hw,color=axes_color)


#initialize waves
dx = 2e-3
waves = []
for n in range(Nwaves):
    waves.append(0)
    
c_pos = []
pts_base = []

A0s = [1 for x in range(Nwaves)]
T0s = [1 for x in range(Nwaves)]
shift0s = [0,pi/2,pi,1.5*pi,pi/4,3*pi/4,5*pi/4,7*pi/4]    #

As = A0s.copy()
Ts = T0s.copy()
shifts = shift0s.copy()

def change_wave(i,A,T,shift):
    #assumes c_pos, pts_base, wave already exist
    global c_pos, waves
    x = 0
    n = 0
    while x < x_max:
        waves[i].modify(n,y=A*sin(2*pi*x/T + shift))
        
        n += 1
        x = n*dx    #more stable way to add

for i,A in enumerate(As):
    x = 0
    n = 0
    if i == 0:
        while x < x_max:
            c_pos.append(vec(x,A*sin(2*pi*x/Ts[i] + shifts[i]),0))
            pts_base.append(vec(x,0,0))
            n += 1
            x = n*dx    #more stable way to add
    else:
        x = 0
        n = 0
        while x < x_max:
            c_pos[n] = vec(x,A*sin(2*pi*x/Ts[i] + shifts[i]), 0)
            
            n += 1
            x = n*dx    #more stable way to add
            
    waves[i] = curve(c_pos)
    waves[i].radius = wave_radfac*x_max
    waves[i].color = wave_colors[i]
    
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
wave0_txt = 'wave1\nA='+str(As[0])+'\n'+'T='+str(Ts[0])+'\n'+'s='+str(shifts[0])
wave0_lbl = label(text=wave0_txt, pos=vec(0,y_max,0),box=False,line=False,
                opacity=0,height=font_h,
                xoffset=-5,align='right',font='monospace',color=wave_colors[0])


eqn_lbl = label(text='y=Asin(2pi x/T + s)',pos=vec(0,y_max/5,0),box=False,line=False,
                opacity=0,height=font_h,
                xoffset=-5,align='right',font='monospace',color=color.white)
sumwave_txt = 'ytot=y1 + y2\n'
sumwave_lbl = label(text=sumwave_txt,pos=vec(0,-y_max/5,0),box=False,line=False,
                height=font_h,opacity=0,
                xoffset=-5,align='right',font='monospace',color=sumwave_color)

wave1_txt = 'wave2\nA='+'{0:.2f}'.format(As[1])+'\n'+'T='+'{0:.2f}'.format(Ts[1])+\
            '\n'+'s='+'{0:.2f}'.format(shifts[1])
wave1_lbl = label(text=wave1_txt, pos=vec(0,-y_max,0),box=False,line=False,
                opacity=0,height=font_h,
                xoffset=-5,align='right',font='monospace',color=wave_colors[1])

#buttons for wave visibility
def toggle_wave_visibility(the_radio):    
    global waves
    
    if '1' in the_radio.text:
        i = 0
    else:
        i = 1
    waves[i].visible = not waves[i].visible
    if waves[i].visible:
        the_radio.text = 'Hide Wave '+str(i+1)
    else:
        the_radio.text = 'Show Wave '+str(i+1)
    

wave1_vis = button(text='Hide Wave 1', pos=scn.title_anchor, 
            bind=toggle_wave_visibility)

wave2_vis = button(text='Hide Wave 2', pos=scn.title_anchor, 
            bind=toggle_wave_visibility)
            

#add sliders for second wave
scn.caption = "Vary 2nd wave's Amplitude: \t\t\t Vary 2nd wave's Period: \t\t\t\t\t Shift the 2nd wave:\n"
def set_A1(s):
    global As, wave1_lbl
    i = 1
    As[i] = s.value
    wave1_txt = 'wave2\nA='+'{0:.2f}'.format(As[i])+'\n'+'T='+'{0:.2f}'.format(Ts[i])+\
                '\n'+'s='+'{0:.2f}'.format(shifts[i])
    wave1_lbl.text=wave1_txt
    change_wave(i,As[i],Ts[i],shifts[i])
    change_sumwave()

A1slider = slider(min=0,value=A0s[1],max=A0s[1]*slide_fac,
                step=slide_step,length=0.333*cw,pos=scn.caption_anchor,bind=set_A1)

scn.append_to_caption(' ')

def set_T1(s):
    global Ts, wave1_lbl
    i = 1
    Ts[i] = s.value
    wave1_txt = 'wave2\nA='+'{0:.2f}'.format(As[i])+'\n'+'T='+'{0:.2f}'.format(Ts[i])+\
                '\n'+'s='+'{0:.2f}'.format(shifts[i])
    wave1_lbl.text=wave1_txt
    change_wave(i,As[i],Ts[i],shifts[i])
    change_sumwave()

T1slider = slider(min=T0s[1]/(2*slide_fac),value=T0s[1],max=T0s[1]*slide_fac,
                step=slide_step,length=0.333*cw,pos=scn.caption_anchor,bind=set_T1)

scn.append_to_caption(' ')

def set_shift1(s):
    global shifts, wave1_lbl
    i = 1
    shifts[i] = s.value
    wave1_txt = 'wave2\nA='+'{0:.2f}'.format(As[i])+'\n'+'T='+'{0:.2f}'.format(Ts[i])+\
                '\n'+'s='+'{0:.2f}'.format(shifts[i])
    wave1_lbl.text=wave1_txt
    change_wave(i,As[i],Ts[i],shifts[i])
    change_sumwave()

T1slider = slider(min=0,value=shift0s[1],max=2*pi,
                pos=scn.caption_anchor,length=0.333*cw,bind=set_shift1)

