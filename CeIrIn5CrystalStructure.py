GlowScript 2.1 VPython
"""
This creates the CeIrInt crystal structure, the heavy-fermion 
superconductor whose superconducting properties I 
studied under pressure from 2007 to 2009.
"""

includeBars = True  #Show bars connecting

a = 4.6662e-10       #m
b = 4.6662e-10       #m
c = 7.5168e-10      #m

#Indium atoms at sites:
Inloc0 = vec(0.5*a,0.5*b,0)
Inloc1 = vec(0,0.5*b,0.3035*c)
Inloc2 = vec(0.5*a,0,0.3035*c)

scene=display(forward=vec(-1,-1,0),up=vec(0,0,1),center=0.5*vec(a,b,c))

#Need sizes
f = 0.375           #scale factor for appearance
rCe = f*2e-10
rIr = f*1e-10
rIn = f*1.5e-10
Cecolor = vec(0.3,0.3,1.0)
Ircolor = color.red
Incolor = color.green
print('Ce: blue, Ir: red, In: green')
print('Tetragonal structure: a = b =',a*1e10,'Angstroms, c =',c*1e10,'Angstroms')

#cycle through four corners
Ces = []
Irs = []
Ins = []

if includeBars:
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
    
