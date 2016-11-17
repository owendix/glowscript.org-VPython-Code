GlowScript 2.2 VPython
#given the probability of success, input by user, visualize how often 
#one of two outcomes occurs (binomial)
#success is visualized as a dart hitting one or another spot on a square field 
#probability of each side is proportional to the radius

#Fivethirtyeight: clinton = 71.4%, trump 28.6%
p = 0.714       #probability of 'success' - the bullseye
p0 = p

left_color = vec(0.3,0.3,1)
right_color = vec(0.85,0,0)       #allow to switch later? may not work with graph
dart_colors = [color.cyan,vec(1,0.2,0.8)]
pmf_color0 = color.green
dart_r = 0.01      #radius of dart
font_h = 24

#create two disks
#success = blue (default, allow to switch with radio)
#always want overall unit radius, area = pi --> pi*r**2 = pi*p: r = sqrt(p)

cw = 1300
ch = 0.4*cw
#will need to adjust center, range for input and description below
scn = canvas(width=cw,height=ch,center=vec(1.7,0.4,0),forward=vec(0,0,-1),
        range=0.75,autoscale=False,userzoom=False,userspin=False)

#board: space of all possible values
len_rl = 0.01   #depth in z direction
ht_rl = 1       #must be for area prob to sum to 1
z_rl=-len_rl/2
right = box(pos=vec(0.5,0.5,z_rl),axis=vec(0,0,-1),length=len_rl,width=ht_rl,
            height=ht_rl,color=right_color,shininess=0)
left = box(pos=vec(p/2,0.5,0.001*dart_r + z_rl),axis=vec(0,0,-1),length=len_rl,
            width=p,height=ht_rl,up=vec(0,1,0),color=left_color,shininess=0)

#labels below, input below
p_label = label(text='{0:0.1f}'.format(p*100)+'%',pos=vec(p/2,-0.05,0),
            align='center',height=font_h,box=False,line=False,opacity=0,
            color=left_color)
pc_label = label(text='{0:0.1f}'.format((1-p)*100)+'%',pos=vec((1+p)/2,-0.05,0),
            align='center',height=font_h,box=False,line=False,opacity=0,
            color=right_color)
instr = label(text='Push spacebar to throw random dart\nReload page to restart',
            pos=vec(0.5,-0.125,0),height=font_h,align='center',box=False,
            line=False,opacity=0,color=color.white)
#make graph
bar_w = 0.3
ax_r = 0.5*dart_r
g_ogn = vec(1.2,0,0)
g_yax = cylinder(pos=g_ogn,axis=vec(0,1,0),radius=ax_r,color=color.white)
g_xax = cylinder(pos=g_ogn,axis=vec(1,0,0),radius=ax_r,color=color.white)
Nl = 0
Nr = 0
Nl0 = 0
barplot_l = box(pos=g_ogn+vec(bar_w,Nl/2,0),axis=vec(0,Nl,0),width=bar_w,
                height=0.1*bar_w,up=vec(0,0,1),color=left_color)
barplot_r = box(pos=g_ogn+vec(1-bar_w,Nr/2,0),axis=vec(0,Nr,0),width=bar_w,
                height=0.1*bar_w,up=vec(0,0,1),color=right_color)   
barplot_l_label = label(text=str(Nl),pos=barplot_l.pos+barplot_l.axis,height=font_h,
                    box=False,line=False,opacity=0,yoffset=0.5*font_h,
                    align='center',color=dart_colors[0])
barplot_r_label = label(text=str(Nr),pos=barplot_r.pos+barplot_r.axis,height=font_h,
                    box=False,line=False,opacity=0,yoffset=0.5*font_h,
                    align='center',color=dart_colors[1])
barplot_yscale = 10

barplot_lp = label(text=0,pos=barplot_l.pos+vec(0,-0.05,0),height=font_h,
                    box=False,line=False,opacity=0,
                    align='center',color=dart_colors[0])
barplot_rp = label(text=0,pos=barplot_r.pos+vec(0,-0.05,0),height=font_h,
                    box=False,line=False,opacity=0,
                    align='center',color=dart_colors[1])

def update_barplot():
    global barplot_l, barplot_r, barplot_yscale, barplot_l_label,barplot_r_label,barplot_lp,barplot_rp

    if Nl >= barplot_yscale or Nr >= barplot_yscale:
        barplot_yscale *= 3
    barplot_l.axis = vec(0,Nl/barplot_yscale,0)
    barplot_l.pos.y = Nl/(2*barplot_yscale)
    barplot_r.axis = vec(0,Nr/barplot_yscale,0)
    barplot_r.pos.y = Nr/(2*barplot_yscale)
    barplot_l_label.text=str(Nl)
    barplot_r_label.text=str(Nr)
    barplot_l_label.pos.y = barplot_l.axis.y
    barplot_r_label.pos.y = barplot_r.axis.y
    lp = 100*Nl/(Nl+Nr)
    barplot_lp.text='{0:0.1f}'.format(lp)+'%'
    barplot_rp.text='{0:0.1f}'.format(100-lp)+'%'
    

def binom_pmf(k,n):
    #k,n integers: 0 <= k <= n
    #0<= p <= 1
    #this one is pretty fast and matches outputs for scipy.stats.binom.pmf
    #uses global variable p
    if p == 0:
        if k == 0:
            return 1
        else:
            return 0
    elif p == 1:
        if k == n:
            return 1
        else:
            return 0
    else:
        if k == 0:
            return (1-p)**n
        elif k == n:
            return p**n
        elif k >= n/2:
            #cancel numerator with x! in denom
            #n!/x! = n*(n-1)*...*(x+1)
            d1 = k
            d2 = n - k
            p1 = p
            p2 = 1-p
        else:
            d1 = n - k
            d2 = k
            p1 = 1-p
            p2 = p
        
        num = n*p2*(p1**d1)

        for q in range(d1+1,n):
            num *= q*p2
            
        for q in range(2,d2+1):
            num /= q
            
        return num

#graph labels for left and right: P(left) = p, P(right) = 1-p
g2_ogn = vec(2.4,0,0)
g2_yax = cylinder(pos=g2_ogn,axis=vec(0,1,0),radius=ax_r,color=color.white)
g2_xax = cylinder(pos=g2_ogn,axis=vec(1,0,0),radius=ax_r,color=color.white)
g2_xax_label = label(text='n_left\nn_left0 = n_left when prob. is adjusted',
                pos=g2_ogn+vec(0.5,-0.125,0),align='center',height=font_h,
                color=left_color,box=False,line=False,opacity=0)
g2_yax_label = label(text='p(n_left | N, p='+'{:0.1f}'.format(p*100)+'%, n_left0=0)',
                pos=g2_ogn+vec(0.5,1.075,0),align='center',height=font_h,
                color=left_color,box=False,line=False,opacity=0)

g2_yticks = [cylinder(pos=g2_ogn+vec(-5*ax_r,0.5,0),axis=vec(10*ax_r,0,0),radius=ax_r,
                color=pmf_color0),
            cylinder(pos=g2_ogn+vec(-5*ax_r,1,0),axis=vec(10*ax_r,0,0),radius=ax_r,
                color=pmf_color0)]
g2_yticklabels = [label(text='0.50',pos=g2_ogn+vec(-0.025,0.5,0),align='right',
                height=font_h,color=pmf_color0,box=False,line=False,opacity=0),
                label(text='1.00',pos=g2_ogn+vec(-0.025,1,0),align='right',
                height=font_h,color=pmf_color0,box=False,line=False,opacity=0)]
g2_xticklabels = [label(text=0,pos=g2_ogn+vec(0,-0.05,0),align='center',height=font_h,
                    color=pmf_color0,box=False,line=False,opacity=0),
                    label(text=0,pos=g2_ogn+vec(0,-0.05,0),align='center',height=font_h,
                    color=left_color,box=False,line=False,opacity=0,visible=False),
                    label(text=0,pos=g2_ogn+vec(0,-0.05,0),align='center',height=font_h,
                    color=pmf_color0,box=False,line=False,opacity=0,visible=False)]

overflow_warning = label(text='Binom_pmf Overflow Error\nAdjust Prob. or Reload Page',
                    pos=g2_ogn + vec(0.5,0.5,0),align='center',height=font_h,
                    color=color.red,box=False,line=False,opacity=0,visible=False)

pmf_pts = []
#as n grows, append to pmf_pts0, if it gets reset, make invisible but don't recreate
pmf = []
Nreset = 0
pmf_h = 0.1*bar_w
pmf_yscale = 1
overflow = False

def set_pmf(f):
    #use current n value, not future
    #future n value probability will be conditional and, bc of independence,
    #will yield just p as the result: bernoulli distribution
    global pmf_pts, pmf, g2_xticklabels, g2_yticklabels, pmf_yscale, overflow_warning
    #maybe have a reset variable, n = 0
    #I'll want only new points, not Nr, Nl
    len_pmf = len(pmf)
    N = Nr + Nl - Nreset
    if N == 0:
        pmf_w = 0
    else:
        pmf_w = 1/(N+1)
    
    #adjust x tick labels: 0 = 0, 1 = Nl, 2 = N
    if N == 1:
        g2_xticklabels[1].visible = True
        g2_xticklabels[2].visible = True
    
    if Nl - Nl0 == 0:
        g2_xticklabels[0].visible = False
    else:
        g2_xticklabels[0].visible = True
        if Nl - Nl0 == N:
            g2_xticklabels[2].visible = False
        else:
            g2_xticklabels[2].visible = True
    
    g2_xticklabels[0].pos.x = g2_ogn.x + pmf_w/2
    g2_xticklabels[1].text = Nl
    g2_xticklabels[1].pos.x = g2_ogn.x + (Nl-Nl0+0.5)*pmf_w
    g2_xticklabels[2].text = N + Nl0
    g2_xticklabels[2].pos.x = g2_ogn.x + 1 - pmf_w/2
    
    
    #adjust yticklabel (just max)
    ysfac = 2
    if N > 5:
        pmf_ymax = f(round(N*p),N)
    else:
        pmf_ymax = 1
    if pmf_ymax < pmf_yscale/ysfac:
        pmf_yscale /= ysfac
    g2_yticklabels[0].text = '{:.2f}'.format(pmf_yscale/2)
    g2_yticklabels[1].text = '{:.2f}'.format(pmf_yscale)
    
    for k in range(N+1):
        if k == Nl-Nl0:
            pmf_color = left_color
        else:
            pmf_color = pmf_color0
        y = f(k,N)
        if y == float('Infinity') or y==float('Inf') or y == float('NaN'):
            for q in pmf:
                q.visible=False
            overflow_warning.visible=True
            break
        y /= pmf_yscale
        if k >= len_pmf:
            pmf_pts.append(vec(k*pmf_w,y,0))
            pmf.append(box(pos=g2_ogn+vec(pmf_pts[k].x+pmf_w/2,pmf_pts[k].y/2,0),
                            axis=vec(0,pmf_pts[k].y,0),up=vec(0,0,1),height=pmf_h,
                            width=pmf_w,color=pmf_color))
        else:
            pmf_pts[k] = vec(k*pmf_w,y,0)
            pmf[k].pos = g2_ogn+vec(pmf_pts[k].x+pmf_w/2,pmf_pts[k].y/2,0)
            pmf[k].axis = vec(0,pmf_pts[k].y,0)
            pmf[k].width = pmf_w
            pmf[k].color = pmf_color
            pmf[k].visible = True    


def reset_pmf():
    global pmf, g2_yticklabels, g2_xticklabels, pmf_yscale, Nl0, Nreset
    global g2_yax_label, overflow, overflow_warning
    Nreset = Nr + Nl
    Nl0 = Nl
    
    if len(pmf) > 0:
        for k in pmf:
            k.visible=False
            
    overflow = False
    overflow_warning.visible = False
    
    pmf_yscale=1
    g2_yax_label.text=text='p(n_left | N,  p='+'{:0.1f}'.format(p*100)+'%, n_left0='+str(Nl)+')'
    g2_yticklabels[0].text='0.50'
    g2_yticklabels[1].text='1.00'
    for i, l in enumerate(g2_xticklabels):
        if i != 0:
            l.visible=False
        l.text=Nl0
        l.pos.x=g2_ogn.x
        
            
#throw dart (small sphere)
darts=[]

def throw_dart(evt):
    global Nl, Nr, darts
    s = evt.key
    if s == ' ':
        v = vec.random()
        r = abs(v.x)
        if r >= p:
            Nr += 1
            dart_color = dart_colors[1]
        else:
            Nl += 1
            dart_color = dart_colors[0]
        theta = pi*v.y
        darts.append(sphere(pos=vec(r,abs(v.y),0),radius=dart_r,
                    color=dart_color))
        update_barplot()
        if not overflow:
            set_pmf(binom_pmf)
    
scn.bind('keydown', throw_dart)


def change_p(s):
    global p, left, p_label, pc_label
    p = s.value
    reset_pmf()

    left.pos.x = p/2
    left.width = p
    p_label.text='{0:0.1f}'.format(p*100)+'%'
    p_label.pos.x = p/2
    pc_label.text='{0:0.1f}'.format((1-p)*100)+'%'
    pc_label.pos.x = (1+p)/2

scn.caption='Adjust probabilities:\n'
slider(min=0.005,value=p0,max=0.995,step=0.001,length=0.36*cw,bind=change_p)

