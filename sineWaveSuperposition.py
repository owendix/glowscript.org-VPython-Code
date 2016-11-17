GlowScript 2.2 VPython
#shows the superposition of two sine waves, one with adjustable amplitude, 
#period, and phase shift

Nwaves = 2      #can include superposition of many waves
wave_radfac = 1/800
sumwave_radfac = 3
sumwave_color = color.white
slide_fac = 4
slide_step = 2e-4
font_h = 18

y_max = 1
x_max = 4


cw = 1100
ch = 0.65*cw
scn = display(width=cw,height=ch,center=vec(0.5*x_max,0,0),range=0.65*x_max,
            userzoom=False,userspin=False)
scn.autoscale=False

#create axes
axes_color = color.gray(0.75)
axes_sw = 1.5*wave_radfac*x_max
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
waves = [0,0]
A0s = [1,1]
T0s = [1,1]
shift0s = [0,pi/2]
wave_colors0 = [color.cyan, color.yellow, color.magenta, color.green, color.red, color.orange, color.blue]
wave_colors = wave_colors0[0:2].copy()
for n in range(2,Nwaves):
    waves.append(0)
    A0s.append(0)   #make the amplitudes 0 after first 2 waves
    T0s.append(1)
    shift0s.append(2*pi*abs(vec.random().x))
    wave_colors.append(wave_colors0[n%len(wave_colors0)])

#A0s = [1 for x in range(Nwaves)]
#T0s = [1 for x in range(Nwaves)]
#shift0s = [0,pi/2,0.6*pi,1.5*pi,pi/4,3*pi/4,5*pi/4,7*pi/4]    #
    
c_pos = []
pts_base = []

As = A0s.copy()
Ts = T0s.copy()
shifts = shift0s.copy()


def change_wave(i,A,T,shift):
    #assumes c_pos, pts_base, wave already exist
    global c_pos, waves
    x = 0
    n = 0
    while x < x_max:
        waves[i].modify(n,y=A*sin(2*pi*x/T - shift))
        
        n += 1
        x = n*dx    #more stable way to add

for i,A in enumerate(As):
    x = 0
    n = 0
    if i == 0:
        while x < x_max:
            c_pos.append(vec(x,A*sin(2*pi*x/Ts[i] - shifts[i]),0))
            pts_base.append(vec(x,0,0))
            n += 1
            x = n*dx    #more stable way to add
    else:
        x = 0
        n = 0
        while x < x_max:
            c_pos[n] = vec(x,A*sin(2*pi*x/Ts[i] - shifts[i]), 0)
            
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
wave0_txt = 'wave1\nA='+str(As[0])+'\n'+'T='+str(Ts[0])+'\n'+'s='+str(shifts[0])
wave0_lbl = label(text=wave0_txt, pos=vec(0,y_max,0),box=False,line=False,
                opacity=0,height=font_h,
                xoffset=-5,align='right',font='monospace',color=wave_colors[0])


eqn_lbl = label(text='y=Asin(2pi x/T - s)',pos=vec(0,y_max/5,0),box=False,line=False,
                opacity=0,height=font_h,
                xoffset=-5,align='right',font='monospace',color=color.white)
if Nwaves <= 2:
    sumwave_txt = 'ytot=y1 + y2\n'
else:
    sumwave_txt =  'ytot=y1 + y2 + ...\n'
sumwave_lbl = label(text=sumwave_txt,pos=vec(0,-y_max/5,0),box=False,line=False,
                height=font_h,opacity=0,
                xoffset=-5,align='right',font='monospace',color=sumwave_color)

wave1_txt = 'wave2\nA='+'{:.2f}'.format(As[1])+'\n'+'T='+'{:.2f}'.format(Ts[1])+'\n'+'s='+'{:.2f}'.format(shifts[1])
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
        wave1_txt = 'wave2\nA='+'{0:.2f}'.format(As[i])+'\n'+'T='+'{0:.2f}'.format(Ts[i])+\
                '\n'+'s='+'{0:.2f}'.format(shifts[i])
        wave1_lbl.text=wave1_txt
    
    change_wave(i,As[i],Ts[i],shifts[i])
    change_sumwave()
    
def set_T(s):
    global Ts, wave1_lbl
    i = s.wavenum
    Ts[i] = s.value
    if i == 1:
        wave1_txt = 'wave2\nA='+'{0:.2f}'.format(As[i])+'\n'+'T='+'{0:.2f}'.format(Ts[i])+\
                '\n'+'s='+'{0:.2f}'.format(shifts[i])
        wave1_lbl.text=wave1_txt
        
    change_wave(i,As[i],Ts[i],shifts[i])
    change_sumwave()
    
def set_shift(s):
    global shifts, wave1_lbl
    i = s.wavenum
    shifts[i] = s.value
    if i == 1:
        wave1_txt = 'wave2\nA='+'{0:.2f}'.format(As[i])+'\n'+'T='+'{0:.2f}'.format(Ts[i])+\
                '\n'+'s='+'{0:.2f}'.format(shifts[i])
        wave1_lbl.text=wave1_txt
        
    change_wave(i,As[i],Ts[i],shifts[i])
    change_sumwave()

#add sliders for second wave and beyond, first wave unshifted
A_sldrs = []
T_sldrs = []
shift_sldrs = []
for i in range(1,Nwaves):
    if i == 1:
        scn.caption = 'Vary wave'+str(i+1)+' Amplitude: \t\t\t\t\t Vary wave'+str(i+1)+' Period: \t\t\t\t\t\t Shift wave'+str(i+1)+':\n'
    else:
        scn.append_to_caption('Vary wave'+str(i+1)+' Amplitude: \t\t\t\t\t Vary wave'+str(i+1)+' Period: \t\t\t\t\t\t Shift wave'+str(i+1)+':\n')
    
    #use A0s[0] for baseline min and max: A0s[2 or higher] start at 0
    A_sldrs.append(slider(min=0,value=A0s[i],max=A0s[0]*slide_fac,
                    pos=scn.caption_anchor,step=slide_step,length=0.333*cw,bind=set_A))
    A_sldrs[i-1].wavenum = i
    #scn.append_to_caption(' ') #doesn't do anything
    
    #use T0s[0] for baseline min and max
    T_sldrs.append(slider(min=T0s[0]/(10*slide_fac),value=T0s[i],max=T0s[0]*slide_fac,
                    pos=scn.caption_anchor,step=slide_step,length=0.333*cw,bind=set_T))
    T_sldrs[i-1].wavenum = i

    #scn.append_to_caption(' ') #doesn't do anything
    
    shift_sldrs.append(slider(min=0,value=shift0s[i],max=2*pi,
                    pos=scn.caption_anchor,length=0.333*cw,bind=set_shift))
    shift_sldrs[i-1].wavenum = i

    scn.append_to_caption('\n')
