GlowScript 2.1 VPython
print('''You win when you reach 1024. Keys: Up, Down, Left, Right
This game is modeled after Gabriele Cirulli's game 2048.
It was made by Owen Dix in October 2016, because the 
original is such a great game, without use of other source 
code, entirely for the fun of building it.''')
#though this uses the wavelength_to_rgb function mentioned below

last_shift = ['right','left','up','down'][floor(3.9999*abs(vec.random().z))]

#create grid, all positive, x and y, z = 0
#each square size = 1-border_thick
#center of first grid point, lower left, is at vec(0,0,0)
the_rate = 875
speed = 1
dt = 2**-6      #keep power of 2 for perfect storage
bdr = 0.2
bdr_color = color.white
bdr_r = bdr/2
bg_color = color.gray(0.6)
font_h = 24 #good for all values
canv_w = 500
canv_h = canv_w
scn_range = 3
#s_idx = 0       #maintain as the lowest index invisible square
squares=[[0,0,0,0],   
        [0,0,0,0],
        [0,0,0,0],
        [0,0,0,0]]
        #[i][j] will be at location pos.x = i, pos.y = j, except when shifting
        #keep others invisible
ps = [13/16, 5/32, 1/32] #probability of getting [2^0, 2^1, 2^2, 2^3] on new gen
#must sum to 1: that's why I do integer multiples of powers of 2-> no comp error
lose = False
won = 0
N = 0

def wavelength_to_rgb(wavelength):
    
    '''This converts a given wavelength of light to an 
    approximate RGB color value, with a false color spectrum
    This scale the wavelength to the visible light spectrum
    
    Based on code by http://www.noah.org/wiki/Wavelength_to_RGB_in_Python,
    which was based on code by Dan Bruton
    http://www.physics.sfasu.edu/astro/color/spectra.html
    '''
    #Scale it so it doesn't have to be in the visible spectrum (false color)
    
    #Convert meters to nanometers
    #wavelength = wavelength*1e9    #already in nanometers
    #lam_r = 750
    #lam_b = 380
    gamma=0.8
    
    if wavelength >= 380 and wavelength <= 440:
        attenuation = 0.3 + 0.7 * (wavelength - 380) / (440 - 380)
        R = ((-(wavelength - 440) / (440 - 380)) * attenuation) ** gamma
        G = 0.0
        B = (1.0 * attenuation) ** gamma
    elif wavelength >= 440 and wavelength <= 490:
        R = 0.0
        G = ((wavelength - 440) / (490 - 440)) ** gamma
        B = 1.0
    elif wavelength >= 490 and wavelength <= 510:
        R = 0.0
        G = 1.0
        B = (-(wavelength - 510) / (510 - 490)) ** gamma
    elif wavelength >= 510 and wavelength <= 580:
        R = ((wavelength - 510) / (580 - 510)) ** gamma
        G = 1.0
        B = 0.0
    elif wavelength >= 580 and wavelength <= 645:
        R = 1.0
        G = (-(wavelength - 645) / (645 - 580)) ** gamma
        B = 0.0
    elif wavelength >= 645 and wavelength <= 750:
        attenuation = 0.3 + 0.7 * (750 - wavelength) / (750 - 645)
        R = (1.0 * attenuation) ** gamma
        G = 0.0
        B = 0.0
    else:
        R = 0.0
        G = 0.0
        B = 0.0
    
    return R, G, B

square_colors=[]    #map wavelength to
nclr = 9
d_color=(750-380)/nclr
for s in [vec(0.333,0.333,0.333),vec(0,0,0)]:
    for n in range(0,nclr+1):
        R,G,B=wavelength_to_rgb(int(n*d_color+380))
        square_colors.append(vec(R,G,B)+s)
        if len(square_colors) == 16:            #num of different possible tiles
            break   #that'll be the second loop
        

scene=display(width=canv_w,height=canv_h,center=vec(1.5,1.5,0),dir=vec(0,0,-1),
        range=scn_range,autoscale=False,userzoom=False,userspin=False,up=vec(0,1,0),
        background=bg_color)

#create grid
grid = []
#place at z = -bdr_r
bdrz = -bdr
s = 0.05*bdr
grid.append(curve(pos=[vec(-0.5-s,-0.5-s,bdrz),vec(3.5+s,-0.5-s,bdrz),
            vec(3.5+s,3.5+s,bdrz),vec(-0.5-s,3.5+s,bdrz),vec(-0.5-s,-0.5-s,bdrz)],
            radius=bdr_r+s,color=bdr_color,shininess=0))
xyvals = [0.5,1.5,2.5]
for X in xyvals:
    grid.append(cylinder(pos=vec(X,-0.5-s,bdrz),axis=vec(0,4+2*s,0),radius=bdr_r,
            color=bdr_color,shininess=0))
for Y in xyvals:
    grid.append(cylinder(pos=vec(-0.5-s,Y,bdrz),axis=vec(4+2*s,0,0),radius=bdr_r,
            color=bdr_color,shininess=0))
        
    
def new_loc():
    
    w = vec.random()        #-1 to 1
    v = floor(3.99999*abs(w.x))
    v += 12
    #start at opposite end, do cumulative distribution
    #use ps
    step = 2*round(abs(w.y)) - 1 #+1, -1
    cd = 0      #cumulative distribution function
    #ps sum to 1
    qs = [1/2,1/4,1/4]  #I want zero chance of being put in last spot
    ij = len(qs)-1
    for i,p in enumerate(qs[:-1]):
        cd += p
        if w.z < cd:
             ij = i
             break
    if last_shift == 'up':
        for j in range(ij,-1,-1):   #there will be a spot toward 0
            for n in range(v,v + step*4,step):
                i = n%4
                if not squares[i][j].visible:
                    return vec(i,j,0)
    elif last_shift == 'down':
        ij = 3 - ij                 #want high prob near 3
        for j in range(ij,4):       #there will be a spot toward j=3
            for n in range(v,v+step*4,step):
                i = n%4
                if not squares[i][j].visible:
                    return vec(i,j,0)
    elif last_shift == 'right':     
        for i in range(ij,-1,-1):   #there will be a spot toward 0
            for n in range(v,v+step*4,step):
                j = n%4
                if not squares[i][j].visible:
                    return vec(i,j,0)
    elif last_shift == 'left':       #'left'
        ij = 3 - ij                 #want high prob near 3
        for i in range(ij,4):       #there will be a spot toward i=3
            for n in range(v,v+step*4,step):
                j = n%4
                if not squares[i][j].visible:
                    return vec(i,j,0)

def new_valpower():
    #get a new value power of 2 proportional to the probability of that value
    v = abs(vec.random().x)
    cd = 0      #cumulative distribution function
    #ps sum to 1
    for i,p in enumerate(ps[:-1]):
        cd += p
        if v < cd:
            return i
    #just in case they don't sum to 1, faster anyway
    return len(ps)-1
        
def make_box():
    #make a box at vector loc, get value from new_valpower()
    #loc.x = {0,1,2,3}
    #loc.y = {0,1,2,3}
    global squares, N, lose
    
    if N < 16:
        loc = new_loc()     #i,j value
        vp = new_valpower()
        #just CHANGE an invisible square, put it in its location, color, label, 
        # and make it visible
        #s_idx is always the lowest index of an invisible square
        squares[loc.x][loc.y].color = square_colors[vp]
        squares[loc.x][loc.y].num.text = str(2**vp)
        squares[loc.x][loc.y].pos=loc
        squares[loc.x][loc.y].stop=loc      #if done here, no need to do it when 
        squares[loc.x][loc.y].add=0         #I do need to maintain stop, add when invisible
        squares[loc.x][loc.y].num.pos=loc
        squares[loc.x][loc.y].visible=True
        squares[loc.x][loc.y].num.visible=True
        squares[loc.x][loc.y].num.exp=vp        #store this for later
        
        N += 1
        if N == 16:
            #check for allowed box combinations
            for i in range(0,4):
                for j in range(0,4):
                    if i < 3:
                        if squares[i][j].num.exp==squares[i+1][j].num.exp:
                            return
                    if i > 0:
                        if squares[i][j].num.exp==squares[i-1][j].num.exp:
                            return
                    if j < 3:
                        if squares[i][j].num.exp==squares[i][j+1].num.exp:
                            return
                    if j > 0:
                        if squares[i][j].num.exp==squares[i][j-1].num.exp:
                            return
            label(text='Game Over\nReload Page To Restart', pos=vec(1.5,1.5,0),
                    height=36, color=color.red, align='center')
            lose = True
    else:
        label(text='Game Over\nReload Page To Restart', pos=vec(1.5,1.5,0),
                height=36, color=color.red, align='center')
        lose = True
       

    
def initialize_squares():
    #make all the boxes, then just make them visible/invisible with time
    for i in range(0,4):
        for j in range(0,4):
            squares[i][j]=box(pos=vec(i,j,0),axis=vec(0,0,1),
                length=bdr_r,width=1-bdr_r,height=1-bdr_r,up=vec(0,1,0),
                color=color.black,shininess=0,visible=False)
            squares[i][j].num = label(text=str(2**0),pos=vec(i,j,0),height=font_h,
                font='monospace',box=False,line=False,opacity=0,
                color=color.white,visible=False)
            squares[i][j].stop = squares[i][j].pos
            squares[i][j].add = 0       #0=no, 1=add to next, 2=add from prev 
            squares[i][j].num.exp=0
            squares[i][j].v = speed
    
    
initialize_squares()    #make the squares
make_box()      #can include up to 4 initial calls to make_box()


n_mv = 0
tmax = 3/speed

new_key=False
shifted = False
processed = True

def keyInput(evt):
    global last_shift, new_key, processed
    s = evt.key
    if processed:
        processed = False
        if s == 'left' or s == 'right' or s == 'up' or s == 'down':        #'left', 'up', 'right','down'
            last_shift = s
            new_key=True
            

scene.bind('keydown', keyInput)


while not lose:
    
    rate(the_rate)
    
    if new_key:
        new_key=False
        #figure out where to shift
        if last_shift == 'right' or last_shift == 'left':
            if last_shift == 'right':
                step=-1
                shift_idx = range(2,-1,step)  #how i will shift
            else:
                step=1
                shift_idx = range(1,4)
            dmax = 0
            for j in range(0,4):
                okto_add = True
                next_sqr = squares[shift_idx[0]-step][j]    #ref to far right sqr
                n_invis = 0     
                if not next_sqr.visible:
                    squares[shift_idx[0]-step][j].stop = vec(shift_idx[-1],j,0)
                    n_invis += 1
                #compare squares to next_sqr, either invisible, same num, or diff num
                for i in shift_idx:
                    if not squares[i][j].visible:
                        squares[i][j].stop = vec(shift_idx[-1] - step*n_invis,j,0)
                        n_invis += 1            
                    else:                   #its a visible square
                        if next_sqr.visible:    #next square visible too
                            if squares[i][j].num.exp == next_sqr.num.exp: #same num
                                if okto_add:
                                    squares[i][j].add = 1
                                    squares[i][j].stop = next_sqr.stop
                                    next_sqr.add = 2
                                    okto_add = False
                                else:
                                    okto_add = True     #reset to true after once
                                    squares[i][j].stop = next_sqr.stop + vec(step,0,0)
                            else:   #visible and different numbers
                                okto_add = True
                                squares[i][j].stop = next_sqr.stop + vec(step,0,0)
                        else:       #its the end of the line
                            squares[i][j].stop = next_sqr.pos
                        d = abs(squares[i][j].stop.x-squares[i][j].pos.x)
                        squares[i][j].v = speed*d   #d is relative to 1
                        #switch order and/or x-->y for other last_shifts
                        if d > dmax:
                            dmax = d
                        next_sqr = squares[i][j]
        
            #determines whether to add a box
            if dmax > 0:
                shifted = True
        
            #now shift
            tmax = 1/speed
            n_mv = 0
            
            #shift and cleanup
            
            if last_shift == 'right':
                #add any boxes
                step = -1
                shift_idx = range(3,-1,step)
            else:       #left
                #for adding boxes
                step = 1
                shift_idx = range(0,4)
                
            
            while n_mv*dt <= tmax:
                rate(the_rate)
                for j in range(0,4):
                    for i in shift_idx[1:]:
                        if squares[i][j].visible:
                            d = -step*(squares[i][j].stop.x - squares[i][j].pos.x)
                            if  d > 0:
                                squares[i][j].pos.x = squares[i][j].pos.x - step*squares[i][j].v*dt
                                squares[i][j].num.pos.x = squares[i][j].pos.x
                n_mv += 1
                
            #add boxes, both left and right
            for j in range(0,4):
                n_invis = 0
                for i in shift_idx:
                    if squares[i][j].add == 1:
                        squares[i][j].add = 0
                        squares[i][j].visible = False
                        squares[i][j].num.visible = False
                        squares[i][j].pos = vec(shift_idx[-1] - step*n_invis,j,0)
                        squares[i][j].stop = squares[i][j].pos
                        squares[i][j].num.pos=squares[i][j].pos
                        n_invis += 1
                        N -= 1
                    elif squares[i][j].add == 2:
                        squares[i][j].add = 0
                        the_exp = squares[i][j].num.exp + 1
                        squares[i][j].color = square_colors[the_exp]
                        squares[i][j].num.exp = the_exp
                        squares[i][j].num.text = str(2**the_exp)
                        if won == 0 and 2**the_exp == 1024:
                            won = 1

                    elif not squares[i][j].visible:
                        squares[i][j].pos = vec(shift_idx[-1] - step*n_invis,j,0)
                        squares[i][j].stop = squares[i][j].pos
                        n_invis += 1

                
            #move pointers
            for j in range(0,4):
                for i in shift_idx:
                    while squares[i][j].pos.x != i or squares[i][j].pos.y != j:
                        new_i = round(squares[i][j].pos.x)
                        new_j = round(squares[i][j].pos.y)
                        tmp_sq = squares[new_i][new_j]  #store anything currently pointed to
                        squares[new_i][new_j] = squares[i][j]
                        squares[i][j] = tmp_sq  


                
        else:   #last_shift either 'up' or 'down'
        
        
            if last_shift == 'up':
                step=-1
                shift_idx = range(2,-1,step)  #how j will change
            elif last_shift == 'down':
                step=1
                shift_idx = range(1,4)  #how j will change
            
            dmax = 0
            for i in range(0,4):
                okto_add = True
                next_sqr = squares[i][shift_idx[0]-step]    #ref to far up sqr
                n_invis = 0     
                if not next_sqr.visible:
                    squares[i][shift_idx[0]-step].stop = vec(i,shift_idx[-1],0)
                    n_invis += 1
                #compare squares to next_sqr, either invisible, same num, or diff num
                for j in shift_idx:
                    if not squares[i][j].visible:
                        squares[i][j].stop = vec(i,shift_idx[-1]-step*n_invis,0)
                        n_invis += 1            
                    else:                   #its a visible square
                        if next_sqr.visible:    #next square visible too
                            if squares[i][j].num.exp == next_sqr.num.exp: #same num
                                if okto_add:
                                    squares[i][j].add = 1
                                    squares[i][j].stop = next_sqr.stop
                                    next_sqr.add = 2
                                    okto_add = False
                                else:
                                    okto_add = True     #reset to true after once
                                    squares[i][j].stop = next_sqr.stop + vec(0,step,0)
                            else:   #visible and different numbers
                                okto_add = True
                                squares[i][j].stop = next_sqr.stop + vec(0,step,0)
                        else:       #its the end of the line
                            squares[i][j].stop = next_sqr.pos
                        d = abs(squares[i][j].stop.y-squares[i][j].pos.y)
                        squares[i][j].v = speed*d   #d is relative to 1
                        #switch order and/or x-->y for other last_shifts
                        if d > dmax:
                            dmax = d
                        next_sqr = squares[i][j]

            
            #determines whether to add a box
            if dmax > 0:
                shifted = True
                
            
            #shift and cleanup
            tmax = 1/speed
            n_mv = 0
            if last_shift == 'up':
                #for adding boxes
                step = -1
                shift_idx = range(3,-1,step)
                
            else:   #last_shift == 'down':
                #for adding boxes
                step = 1
                shift_idx = range(0,4)
            
            #shift
            while n_mv*dt <= tmax:
                rate(the_rate)
                for i in range(0,4):
                    for j in shift_idx[1:]:
                        if squares[i][j].visible:
                            d = -step*(squares[i][j].stop.y - squares[i][j].pos.y)
                            if d > 0:
                                squares[i][j].pos.y = squares[i][j].pos.y - step*squares[i][j].v*dt
                                squares[i][j].num.pos.y = squares[i][j].pos.y
                n_mv += 1
            
            
            #add boxes, both up and down
            for i in range(0,4):
                n_invis = 0
                for j in shift_idx:
                    #do adds, will probably need to make it move right
                    if squares[i][j].add == 1:
                        squares[i][j].add = 0
                        squares[i][j].visible = False
                        squares[i][j].num.visible = False
                        squares[i][j].pos = vec(i,shift_idx[-1] - step*n_invis,0)
                        squares[i][j].stop = squares[i][j].pos
                        squares[i][j].num.pos=squares[i][j].pos
                        n_invis += 1
                        N -= 1
                    elif squares[i][j].add == 2:
                        squares[i][j].add = 0
                        the_exp = squares[i][j].num.exp + 1
                        squares[i][j].color = square_colors[the_exp]
                        squares[i][j].num.exp = the_exp
                        squares[i][j].num.text = str(2**the_exp)
                        if won == 0 and 2**the_exp == 1024:
                            won = 1

                    elif not squares[i][j].visible:
                        squares[i][j].pos = vec(i, shift_idx[-1] - step*n_invis,0)
                        squares[i][j].stop = squares[i][j].pos
                        n_invis += 1
                        
            
            #move pointers
            for i in range(0,4):
                for j in shift_idx:
                    while squares[i][j].pos.x != i or squares[i][j].pos.y != j:
                        new_i = round(squares[i][j].pos.x)
                        new_j = round(squares[i][j].pos.y)
                        tmp_sq = squares[new_i][new_j]  #store anything currently pointed to
                        squares[new_i][new_j] = squares[i][j]
                        squares[i][j] = tmp_sq 
                        
        if won == 1:
            won = 2
            win_label = label(text='You Win\nClick To Continue\nReload To Restart', 
                                pos=vec(1.5,1.5,0),
                                height=36, color=color.red, align='center')
            scene.waitfor('click')
            win_label.visible=False
        
        if shifted:       
            make_box()
            
        shifted = False
        processed = True
        
