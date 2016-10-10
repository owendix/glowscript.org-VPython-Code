GlowScript 2.1 VPython
#plot Earth's g as a function of distance from the center, from the radius out
plot_moon = True

Re = 6.371e6    #m
r_i = Re
r = r_i
r_moon_orbit = [3.626e8,3.85e8,4.054e8] #meters, perigee, mean, apogee
if not plot_moon:
    r_max = 4*Re
else:
    r_max = r_moon_orbit[-1]*1.05

G = 6.674e-11   #Nm^2/kg^2
Me = 5.972e24   #kg
cw = 800
ch = (9/16)*cw
#graph stuff
g_graph = graph(width=cw,height=ch,align='left',xmin=0,ymin=0)
g_series = series(graph=g_graph,label='g (m/s^2)')


dr = (r_max - r)/100
g = G*Me/r**2

while r < r_max:
    g = G*Me/r**2
    g_series.plot(r,g)
    r = r + dr


h_iss = 4e5    #meters elevation of iss orbit
h_atm = 1e5         #meters edge of atmosphere
h_plane = 9144      #meters, plane elevation (30,000 feet)
h_hubble = 5.59e5   #meters, hubble orbit height
hs = [h_plane, h_atm, h_iss, h_hubble]
h_strs = ['sea level', 'plane', 'atm edge', 'iss', 'hubble']
r_keyvals = [Re]
for h in hs:        #construct key values
    r_keyvals.append(Re+h)

the_colors = [color.red,color.orange,color.yellow,color.green,color.blue,color.magenta]
g_gdots = []
for i, h_str in enumerate(h_strs):
    g_gdots.append(gdots(graph=g_graph,label=h_str,color=the_colors[i]))
    r = r_keyvals[i]
    g = G*Me/r**2
    g_gdots[i].plot(r,g)

if plot_moon:
    print('Moon\'s orbit of Earth:')
    r_strs = ['perigee', 'mean orbit', 'apogee']
    g_gdots.append(gdots(graph=g_graph,label='moon\'s orbit', color=the_colors[-1]))
    rg = []
    for i, r in enumerate(r_moon_orbit):
        g = G*Me/r**2
        rg.append([r,g])
        print(r_strs[i],': g =',g, 'm/s^2')
    g_gdots[-1].plot(rg)
