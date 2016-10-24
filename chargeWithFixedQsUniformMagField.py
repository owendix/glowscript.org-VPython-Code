GlowScript 2.1 VPython


the_rate = 50000

#free particle's properties, the system
sys_pos_init = vec(-1.5,1,0)
sys_v_init = vec(0,0,0)
sys_m = 1.5
sys_q = 2
sys_color = color.red
retain_trail = 2000

#create magnetic field: uniform
B = vec(-1,4,-3)        #Tesla, if you wish (which would be very strong)
#fixed charge properties
fixedQs_pos= [vec(0,0,0), vec(-3,2,0), vec(0,-3,4)]
fixedQs_color = [color.magenta, color.yellow, color.cyan]
fixedQs_q = [1,-2,3]


cw = 1000
ch = (9/16)*cw
scn = canvas(width=cw,height=ch,background=color.gray(0.1))

#create fixed charges
rQ = 0.2
fixedQs = []
for i, p in enumerate(fixedQs_pos):
    fixedQs.append(sphere(pos=p,radius=rQ,color=fixedQs_color[i]))
    fixedQs[-1].q = fixedQs_q[i]
    fixedQs[-1].lbl = label(text=str(fixedQs_q[i]),pos=p,xoffset=10,yoffset=10,
                    height=18, box=False,line=False,opacity=0,
                    color=fixedQs_color[i],align='center')

#create system particle
sys = sphere(pos=sys_pos_init, radius=rQ, color=sys_color,make_trail=True,
    retain=retain_trail)
sys.v = sys_v_init
sys.m = sys_m
sys.q = sys_q
sys.lbl = label(text=str(sys.q),pos=sys.pos,xoffset=10,yoffset=10,
                    height=18, box=False,line=False,opacity=0,
                    color=sys_color,align='center')

#k = 8.988e9 #Nm^2/C^2
k = 1

t = 0
dt = 2**-13
fscale = 1
swscale = 0.01
n = 0
n_print = the_rate
while True:
    rate(the_rate)
    sys.pos = sys.pos + sys.v*dt
    sys.lbl.pos = sys.pos
    
    #sum up the forces
    Fnet = sys.q*cross(sys.v,B)
    #if n%n_print == 0:
    #    print('fq.q, F, |F|, Fnet, |Fnet|')
    for fq in fixedQs:
        rvec = sys.pos - fq.pos
        F = (k*sys.q*fq.q/mag2(rvec))*hat(rvec)
        Fnet = Fnet + F
    #    if n%n_print == 0:
    #        print(fq.q, F, mag(F), Fnet, mag(Fnet))
    #if n%n_print == 0: 
    #    scn.waitfor('click')
    
    if t == 0:
        #Fvec = arrow(pos=sys.pos, axis=fscale*Fnet, 
        #    shaftwidth=swscale*fscale*mag(Fnet),
        #    color=sys_color)
        Flabel = label(text='Fnet=('+'{0:.1f}'.format(Fnet.x)+','+
            '{0:.1f}'.format(Fnet.y)+','+'{0:.1f}'.format(Fnet.z)+')',
            pos=scn.center+vec(-0.9*scn.range,-0.9*scn.range,0),
            height=18,box=False,line=False,opacity=0,color=sys_color,
            align='center')
    elif n%n_print == 0:
        #Fvec.pos = sys.pos
        #Fvec.axis=fscale*Fnet
        #Fvec.shaftwidth=swscale*mag(Fvec.axis)
        Flabel.text='Fnet=('+'{0:.1f}'.format(Fnet.x)+','+\
            '{0:.1f}'.format(Fnet.y)+','+'{0:.1f}'.format(Fnet.z)+')'
        Flabel.pos=scn.center+vec(-0.9*scn.range,-0.9*scn.range,0)
    #print(Fnet)
    
    #scn.waitfor('click')
        
    sys.v = sys.v + Fnet*(dt/sys.m)
    n += 1
    t = n*dt
    
