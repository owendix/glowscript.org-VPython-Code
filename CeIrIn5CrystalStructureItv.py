GlowScript 2.1 VPython
"""
This creates the CeIrInt crystal structure, the heavy-fermion 
superconductor whose superconducting properties I 
studied under pressure from 2007 to 2009.
"""

the_rate = 100
dt = 1e-4
ang0 = radians(20)
canv_w = 1200
a = 4.6662e-10       #m
b = 4.6662e-10       #m
c = 7.5168e-10      #m

the_center=0.5*vec(a,b,c)

scene=display(width=canv_w,height=(9/16)*canv_w,
        forward=vec(-1,-1,0),up=vec(0,0,1),center=the_center)

scene.title='CeIrIn_5 Crystal Structure\n\n'

running = True

def Run(b):
    global running
    running = not running
    if running: 
        b.text = "Pause"
    else: 
        b.text = "Run"
    
button(text="Pause", pos=scene.title_anchor, bind=Run)

include_bars = True  #Show bars connecting

#Indium atoms at sites:
Inloc0 = vec(0.5*a,0.5*b,0)
Inloc1 = vec(0,0.5*b,0.3035*c)
Inloc2 = vec(0.5*a,0,0.3035*c)        


#Need sizes
f = 0.375           #scale factor for appearance
rCe = f*2e-10
rIr = f*1e-10
rIn = f*1.5e-10
Cecolor = vec(0.3,0.3,1.0)
Ircolor = color.red
Incolor = color.green


#cycle through four corners
Ces = []
Irs = []
Ins = []


#tetragonal structure
r = 0.02*a
structs=[]
origin=vec(0,0,0)
opac=0.25
structs.append(cylinder(pos=origin,axis=vec(a,0,0),radius=r,color=color.white,
    opacity=opac))
structs.append(cylinder(pos=origin,axis=vec(0,b,0),radius=r,color=color.white,
    opacity=opac))
structs.append(cylinder(pos=origin,axis=vec(0,0,c),radius=r,color=color.white,
    opacity=opac))
bar_shifts_a=[vec(0,b,0),vec(0,0,c),vec(0,b,c)]#,vec(0,0,0.5*c),vec(0,b,0.5*c)]
bar_shifts_b=[vec(a,0,0),vec(0,0,c),vec(a,0,c)]#,vec(0,0,0.5*c),vec(a,0,0.5*c)]
bar_shifts_c=[vec(a,0,0),vec(0,b,0),vec(a,b,0)]
for i in range(len(bar_shifts_a)):
    structs.append(structs[0].clone(pos=bar_shifts_a[i]))
    structs.append(structs[1].clone(pos=bar_shifts_b[i]))
    if i<len(bar_shifts_c):
        structs.append(structs[2].clone(pos=bar_shifts_c[i]))


see_bars = True

def See_Bars(sb_radio):
    global see_bars, structs     #make things global if its changed in here
    see_bars = not see_bars
    for b in structs:
        b.visible = see_bars
    if see_bars:
        sb_radio.text = "Hide Guides"
    else:
        sb_radio.text = "Show Guides"
      
button(text="Hide Guides", pos=scene.title_anchor, bind=See_Bars)

lat_shifts = [vec(0,0,0),vec(a,0,0),vec(0,b,0),vec(a,b,0)]
#Create Ceriums and Iridiums
for s in lat_shifts:
    #Create Ceriums (bottom)
    Ces.append(sphere(pos=s,radius=rCe,color=Cecolor))
    #Create Iridiums (middle)
    Irs.append(sphere(pos=s+vec(0,0,0.5*c),radius=rIr,color=Ircolor))
    #Create Ceriums (top)
    Ces.append(sphere(pos=s+vec(0,0,c),radius=rCe,color=Cecolor))

#Create Indiums
#bottom
Ins.append(sphere(pos=Inloc0,radius=rIn,color=Incolor))
#top
Ins.append(sphere(pos=Inloc0+vec(0,0,c),radius=rIn,color=Incolor))
#lower sides
Ins.append(sphere(pos=Inloc1,radius=rIn,color=Incolor))
Ins.append(sphere(pos=Inloc1+vec(a,0,0),radius=rIn,color=Incolor))
Ins.append(sphere(pos=Inloc2,radius=rIn,color=Incolor))
Ins.append(sphere(pos=Inloc2+vec(0,b,0),radius=rIn,color=Incolor))
#upper sides
Ins.append(sphere(pos=vec(0,b,c)-Inloc1,radius=rIn,color=Incolor))
Ins.append(sphere(pos=vec(a,b,c)-Inloc1,radius=rIn,color=Incolor))
Ins.append(sphere(pos=vec(a,0,c)-Inloc2,radius=rIn,color=Incolor))
Ins.append(sphere(pos=vec(a,b,c)-Inloc2,radius=rIn,color=Incolor))

#shift atoms down
atoms=[]
atoms.extend(Ces)
atoms.extend(Irs)
atoms.extend(Ins)

scene.caption = "Vary the rotation speed: \n"
omega0 = 60
omega = omega0
def set_omega(s):
    global omega
    omega = s.value
    
#default at caption anchor
slider(min=-300, value=omega0, max=300, length=canv_w, bind=set_omega)

#Initial rotation
for s in atoms:
    s.rotate(angle=ang0, axis=-scene.forward, origin=the_center)
for s in structs:
    s.rotate(angle=ang0, axis=-scene.forward, origin=the_center)
scene.range=0.75*c+2*rCe

print('Ce: blue, Ir: red, In: green')
print('Tetragonal structure: a = b =',a*1e10,'Angstroms, c =',c*1e10,'Angstroms')
print('CeIrIn_5 is an unconventional superconductor at T < 0.4K, at 0 pressure and magnetic field')
while True:
    rate(the_rate)
    if running:
        for s in atoms:
            s.rotate(angle=omega*dt, axis=scene.up, origin=the_center)
        for s in structs:
            s.rotate(angle=omega*dt, axis=scene.up, origin=the_center)
