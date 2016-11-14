GlowScript 2.2 VPython
#initial params
x_max = 6.1375
x_min = 0.7
y_max = 2.75
y_min = -2.75
dx = 2**-5
mv_eps = 0.25*dx
crv_r = 2**-6
crv_color = color.orange
font_h = 18
pt1_color = color.cyan
pt2_color = color.yellow
pt_color = vec(0.8,0.3,0.8)
slope_ln_color = color.green
highlight_color=color.white
highlight_opac = 0.25
instr_color = color.white
pt_rfac = 4


#set scene
cw = 900
ch = 0.55*cw
#change center, range
scn = canvas(width=cw, height=ch, center=vec(0.33*x_max,(y_min+y_max)/2,0),
            range=0.55*x_max, forward=vec(0,0,-1), userzoom=False, 
            userspin=False)
scn.autoscale = False
#make axes
xaxes_color = color.white
yaxes_color = color.white
axes_sw = 1.5*crv_r
axes_hl = 5*axes_sw
axes_hw = 0.5*axes_hl
axes_pad = 1.2
ogn = vec(0,0,0)
ax_xp = arrow(pos=ogn,axis=axes_pad*vec(x_max,0,0),shaftwidth=axes_sw,
            headlength=axes_hl,headwidth=axes_hw,color=xaxes_color)
ax_yp = arrow(pos=ogn,axis=axes_pad*vec(0,y_max,0),shaftwidth=axes_sw,
            headlength=axes_hl,headwidth=axes_hw,color=yaxes_color)
ax_yn = arrow(pos=ogn,axis=axes_pad*vec(0,-y_max,0),shaftwidth=axes_sw,
            headlength=axes_hl,headwidth=axes_hw,color=yaxes_color)
y_lbl = label(text='y',pos=ax_yp.pos+ax_yp.axis,xoffset=-5,yoffset=-5,height=font_h,
            box=False,line=False,opacity=0,align='right',color=yaxes_color)
x_lbl = label(text='x',pos=ax_xp.pos+ax_xp.axis,yoffset=-5,height=font_h,
            box=False,line=False,opacity=0,align='center',color=xaxes_color)            
                
#make function
fxs_pos = []
dfxs_pos = []
x = x_min
a = 0.3
b = 3.5
c = -1
d = 4
def f(x):
    return a*(x-b)**3 + c*x + d

def df(x):
    return 3*a*(x-b)**2 + c

while x < x_max:
    fxs_pos.append(vec(x,f(x),0))
    x += dx    #more stable way to add
    

#plot function
fx = curve(fxs_pos)
fx.radius=crv_r
fx.color=crv_color

#give instructions:
instr=label(text='Click for instantaneous slope\nClick and drag for average slope',
            pos=scn.center+vec(-1.2*scn.range,y_min,0),box=False,line=False,
            align='center',opacity=0,color=instr_color,height=font_h)
#point/region markers
highlight = box(pos=vec(0,0,0),axis=vec(0,1.5*(y_max-y_min),0),height=crv_r,
            up=vec(0,0,1),width=dx,color=highlight_color,opacity=highlight_opac,
            visible=False)
pt1 = sphere(pos=vec(0,0,0),radius=pt_rfac*crv_r,color=pt1_color,visible=False)
pt2 = sphere(pos=vec(0,0,0),radius=pt_rfac*crv_r,color=pt2_color,visible=False)
pt = sphere(pos=vec(0,0,0),radius=pt_rfac*crv_r,color=pt_color,visible=False)
#line
slope_ln = cylinder(pos=vec(0,0,0),axis=vec(1,0,0),radius=crv_r,color=slope_ln_color,
                    visible=False)
#labels for pt/line values
pt1_lbl = label(text='p1 = ('+'{:+.2f}'.format(pt1.pos.x)+', '+'{:+.2f}'.format(pt1.pos.y)+', 0)',
            pos=scn.center+vec(-1.65*scn.range,y_max,0),box=False,line=False,
            align='left',opacity=0,color=pt1_color,height=font_h,visible=False)
pt2_lbl = label(text='p2 = ('+'{:+.2f}'.format(pt2.pos.x)+','+'{:+.2f}'.format(pt2.pos.y)+', 0)',
            pos=scn.center+vec(-1.65*scn.range,0.85*y_max,0),box=False,line=False,
            align='left',opacity=0,color=pt2_color,height=font_h,visible=False)
pt_lbl = label(text='p = ('+'{:+.2f}'.format(pt.pos.x)+','+'{:+.2f}'.format(pt.pos.y)+', 0)',
            pos=scn.center+vec(-1.65*scn.range,y_max,0),box=False,line=False,
            align='left',opacity=0,color=pt_color,height=font_h,visible=False)
m_lbl = label(text='m = ('+'{:+.2f}'.format(pt2.pos.y)+' - '+
            '{:+.2f}'.format(pt1.pos.y)+')/('+'{:+.2f}'.format(pt2.pos.x)+' - '+
            '{:+.2f}'.format(pt1.pos.x)+')=',
            pos=scn.center+vec(-1.65*scn.range,0.7*y_max,0),box=False,line=False,
            align='left',opacity=0,color=slope_ln.color,height=font_h,
            visible=False)
ln_lbl = label(text='y = mx + b',pos=scn.center+vec(-1.65*scn.range,0,0),box=False,
            line=False,align='left',opacity=0,color=slope_ln.color,height=font_h,
            visible=False)

#click drag: average at t1, t2
#click no drag (below threshold): instantaneous, at t
mark = False
drag = False
t1 = 0
t2 = 0
#add vertical white lines and semi-opaque boxes later
scn.bind("mousedown", def ():
    global mark, t1, pt, pt1, pt2, highlight, slope_ln, pt1_lbl, pt2_lbl, pt_lbl, m_lbl, ln_lbl
    p = scn.mouse.pos
    pt.visible=False
    pt1.visible=False
    pt2.visible=False
    highlight.visible=False
    slope_ln.visible=False
    pt1_lbl.visible=False
    pt2_lbl.visible=False
    pt_lbl.visible=False
    m_lbl.visible=False
    ln_lbl.visible=False
    if p.x >= fxs_pos[0].x and p.x <= fxs_pos[-1].x:
        mark = True
        t1 = p.x
        highlight.pos.x=t1
        highlight.width=dx
        highlight.visible=True
)

scn.bind("mousemove", def ():
    global drag, highlight
    if mark and abs(scn.mouse.pos.x - t1) >= mv_eps: # mouse button is down
        drag = True
        if scn.mouse.pos.x > fxs_pos[-1].x:
            t_lim = fxs_pos[-1].x
        elif scn.mouse.pos.x < fxs_pos[0].x:
            t_lim = fxs_pos[0].x
        else:
            t_lim = scn.mouse.pos.x
            
        highlight.pos.x = (t_lim + t1)/2
        highlight.width = abs(t_lim - t1)
)

scn.bind("mouseup", def ():
    nonlocal mark, drag, t2, slope_ln, pt, pt1, pt2, highlight, m_lbl,pt_lbl,pt1_lbl,pt2_lbl,ln_lbl
    #mark = False    #really?
    #drag = False
    t2 = scn.mouse.pos.x
    highlight.visible=False
    if not mark:
        return
    else:
        if drag:
            if t2 <= fxs_pos[0].x:
                t2 = fxs_pos[0].x
            elif t2 >= fxs_pos[-1].x:
                t2 = fxs_pos[-1].x
        #be sure it didn't be 'dragged' but outside of the function domain
        if abs(t1 - t2) < mv_eps:
            drag = False
        if not drag:
            #go with t1, instantaneous slope
            y = f(t1)
            m = df(t1)
            b = y - m*t1
            pt.pos=vec(t1,f(t1),0)
            pt.visible = True
            #drag slope_ln, instantaneous
            #change color to pt_color
            slope_ln.pos=vec(0,b,0)
            slope_ln.axis = vec(x_max,m*x_max,0)  #make it span the function
            slope_ln.color=pt_color
            slope_ln.visible=True
            pt_lbl.text='p = ('+'{:+.2f}'.format(pt.pos.x)+', '+'{:+.2f}'.format(pt.pos.y)+', 0)'
            pt_lbl.visible=True
            m_lbl.text=text='m = df(x)/dx = '+'{:+.2f}'.format(m)
            m_lbl.color=slope_ln.color
            m_lbl.visible=True
            ln_lbl.text='y = '+'{:.2f}'.format(m)+'x + '+'{:.2f}'.format(b)
            ln_lbl.color=slope_ln.color
            ln_lbl.visible=True

        else:
            y1 = f(t1)
            y2 = f(t2)
            pt1.pos = vec(t1,y1,0)
            pt2.pos = vec(t2,y2,0)
            pt1.visible = True
            pt2.visible = True
            m = (y2 - y1)/(t2-t1)
            b = y1 - m*t1
            slope_ln.pos=vec(0,b,0)
            slope_ln.axis = vec(x_max,m*x_max,0)
            slope_ln.color = slope_ln_color
            slope_ln.visible = True
            pt1_lbl.text='p1 = ('+'{:+.2f}'.format(pt1.pos.x)+', '+'{:+.2f}'.format(pt1.pos.y)+', 0)'
            pt1_lbl.visible=True
            pt2_lbl.text='p2 = ('+'{:+.2f}'.format(pt2.pos.x)+', '+'{:+.2f}'.format(pt2.pos.y)+', 0)'
            pt2_lbl.visible=True
            m_lbl.text=text='m = ('+'{:+.2f}'.format(pt2.pos.y)+' - '+'{:+.2f}'.format(pt1.pos.y)+')/('+'{:+.2f}'.format(pt2.pos.x)+' - '+'{:+.2f}'.format(pt1.pos.x)+')\nm = '+'{:+.2f}'.format(m)
            m_lbl.color=slope_ln.color
            m_lbl.visible=True
            ln_lbl.text='y = '+'{:.2f}'.format(m)+'x + '+'{:.2f}'.format(b)
            ln_lbl.color=slope_ln.color
            ln_lbl.visible=True

        
        drag = False
        mark = False
        
)

