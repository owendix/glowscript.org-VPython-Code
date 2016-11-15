GlowScript 2.2 VPython
#solves a sudoku, based off of your input values, using a brute force method
font_h = 36
font_w2 = 0.175   #half a character length in real space, not pixels
bbdr_r = 0.05
lbdr_r = 0.5*bbdr_r
bdr_color = color.black
valid_input_num_color = color.blue
logic_num_color = color.green   #not yet used, not really needed: it's fast
num_color = color.black
invalid_color = color.red

#board[i][j] will be the ith row and jth column located at pos=vec(j,i,0)
#0th row at bottom of board when displayed, 0th column at left
board = []
nums = []
row=[0,-1,-2,-3,-4,-5,-6,-7,-8]
board.append(row.copy())
summand = -10
for i in range(1,9):
    for j in range(9):
        row[j] += summand
    board.append(row.copy())
'''Initial values for board
#all non-positive, concat: -row#col#, or -(10*row# + col#)
#this lets all initial entries be unique from each other and future entries
#which lets me use list(set(row)) or list(set(col)) for comparisons
#SET IS NOT DEFINED, NEED TO RETHINK THIS
board = [[0, -1, -2, -3, -4, -5, -6, -7, -8],
[-10, -11, -12, -13, -14, -15, -16, -17, -18],
[-20, -21, -22, -23, -24, -25, -26, -27, -28],
[-30, -31, -32, -33, -34, -35, -36, -37, -38],
[-40, -41, -42, -43, -44, -45, -46, -47, -48],
[-50, -51, -52, -53, -54, -55, -56, -57, -58],
[-60, -61, -62, -63, -64, -65, -66, -67, -68],
[-70, -71, -72, -73, -74, -75, -76, -77, -78],
[-80, -81, -82, -83, -84, -85, -86, -87, -88]
]
'''    
cw = 600
ch = cw
scn = canvas(width=cw,height=cw,center=vec(4,3.65,0),forward=vec(0,0,-1),
        range=5,autoscale=False,userzoom=False,userspin=False,
        background=color.gray(0.99))

#initialize labels, board and labels will be changed
#make sure they stay the same once visible but only board is used for 
#checking validity and solving
nums = [0,0,0,0,0,0,0,0,0]
row=[0,0,0,0,0,0,0,0,0]
for i in range(9):
    for j in range(9):
        row[j] = label(text='',pos=vec(j,i,0),height=font_h,align='center',
                box=False,line=False,opacity=0,font='monospace',
                color=valid_input_num_color,visible=True)
    nums[i] = row.copy()

def init_entry_at(i,j):
    #returns single initial entry consistent with above
    #i = row#, j = col#
    return -(10*i + j)
    
def print_board():
    #print(...,sep='') doesn't work in glowscript
    for i in range(8,-1,-1):
        """
        if i == 0:
            print('[',board[i][0],', ',board[i][1],', ',board[i][2],', ',
            board[i][3],', ',board[i][4],', ',board[i][5],', ',board[i][6],', ',
            board[i][7],', ',str(board[i][8])+']')
        else:
        """
        print(board[i])


#create blank initial game board
#create grid
grid = []
#place at z = 0
bdrz = 0
s = 0
grid.append(curve(pos=[vec(-0.5-s,-0.5-s,bdrz),vec(8.5+s,-0.5-s,bdrz),
            vec(8.5+s,8.5+s,bdrz),vec(-0.5-s,8.5+s,bdrz),vec(-0.5-s,-0.5-s,bdrz)],
            radius=bbdr_r+s,color=bdr_color,shininess=0))
xyvals = range(0.5,9,1)
for X in xyvals:
    if X == 2.5 or X == 5.5:
        bdr_r = bbdr_r
    else:
        bdr_r = lbdr_r
    grid.append(cylinder(pos=vec(X,-0.5-s,bdrz),axis=vec(0,9+2*s,0),radius=bdr_r,
            color=bdr_color,shininess=0))
for Y in xyvals:
    if Y == 2.5 or Y == 5.5:
        bdr_r = bbdr_r
    else:
        bdr_r = lbdr_r
    grid.append(cylinder(pos=vec(-0.5-s,Y,bdrz),axis=vec(9+2*s,0,0),radius=bdr_r,
            color=bdr_color,shininess=0))

#instructions, including reload page to retry

#equipment for checking validity
#use set of row and column
temp_list = [0,0,0,0,0,0,0,0,0]
def get_boxlist(i,j):
    #given an entry, return a single list of all entries in its box
    #starting at top left of box going to bottom right
    #made difficult because slicing doesn't work for 2 indices, with glowscript
    #board[0:3][0:3] returns same as board[:][0:3] and as board[0:3][:]
    global temp_list
    n_i0 = 3*int(i/3)  #start of box, row
    n_j0 = 3*int(j/3)  #start of box, col
    n = 0
    for q in range(n_i0, n_i0+3):
        for r in range(n_j0, n_j0+3):
            temp_list[n] = board[q][r]
            n += 1
            
    return temp_list

def set(a_list):
    #dict() doesn't seem to work
    #set doesn't work so I have to improvise with dict
    b = {}
    for a in a_list:
        b[a] = 0
    return b

def is_valid_at(i,j):
    #checking for repeats in row, col, box for i,j
    #glowscript doesn't allow sets, or valid slicing with 2D arrays
    global temp_list
    if len(set(board[i])) < 9:
        return False
    for q in range(9):
        temp_list[q]=board[q][j]
    if len(set(temp_list)) < 9:
        return False
    if len(set(get_boxlist(i,j))) < 9:
        return False
    
    return True
    
def is_valid_board():
    #checking for repeats
    #in each row, column, and unique box
    global temp_list
    for i_row in range(9):
        if len(set(board[i_row])) < 9:
            #not unique, has repeats
            print('i_row:',i_row,'Row error')
            return False
    #in each column
    for i_col in range(9):
        for i in range(9):
            temp_list[i]=board[i][i_col]
        if len(set(temp_list)) < 9:
            #not unique, has repeats
            print('i_col:',i_col,'Column error')
            return False
    #in each box
    for i_box in range(0,9,3):
        for j_box in range(0,9,3):
            if len(set(get_boxlist(i_box,j_box))) < 9:
                print('i_row,i_col:',i_row,i_col,'Box error')
                return False
    
    return True

#allow for input
#board[row][col], nums pos: vec(col,row,0), keep board and nums text same if visible
active_ij = [-1,-1]
cursor = label(text='|',pos=vec(active_ij[1],active_ij[0],0),height=font_h,
            font='monospace',box=False,line=False,opacity=0,align='center',
            color=valid_input_num_color,visible=False)
invalid_warning = label(text='Invalid input',pos=vec(4,-1,0),height=font_h,
            box=False,line=False,opacity=0,align='center',color=invalid_color,
            visible=False)
ok_to_input = True #turn this off when solver is used, reload to start over

def input_click():
    global active_ij, cursor
    p = scn.mouse.pos
    if ok_to_input:
        if p.x > -0.5 and p.x < 8.5 and p.y > -0.5 and p.y < 8.5:
            j = round(p.x)    #correct, not switched, need global variables for these
            i = round(p.y)
            active_ij = [i,j]
            the_num = nums[active_ij[0]][active_ij[1]]
            cursor.pos = vec(active_ij[1],active_ij[0],0)
            if len(the_num.text) > 0:
                cursor.pos.x += font_w2
            cursor.visible=True
        else:
            active_ij = [-1,-1]
            cursor.visible=False
            #put vertical bar there, don't make it blink (for now), purely aesthetic
            #have to move it if active_ij moves with arrow
    
scn.bind('mouseup', input_click)

def key_input(evt):
    global active_ij, board, cursor, nums, invalid_warning
    s = evt.key
    #how to handle backspace, enter, clicking away
    if active_ij[0] > -1:
        the_num = nums[active_ij[0]][active_ij[1]]
        if s == 'backspace' and len(the_num.text) > 0:
            the_num.text = ''
            board[active_ij[0]][active_ij[1]] = init_entry_at(active_ij[0],active_ij[1])
            cursor.pos.x = int(cursor.pos.x)
            is_val = is_valid_at(active_ij[0],active_ij[1])
            if not is_val:
                the_num.color=invalid_color
                invalid_warning.visible=True
            else:
                the_num.color=valid_input_num_color
                invalid_warning.visible=False
        elif s=='\n' or s=='up' or s=='down' or s=='left' or s=='right' or s=='enter' or s=='return':
            if s == 'down' or s == '\n' or s == 'enter' or s == 'return':
                active_ij[0] = (active_ij[0] + 9 - 1) % 9
            elif s=='right':    #don't allow tab, it activates address bar
                active_ij[1] = (active_ij[1] + 1) % 9
            elif s == 'up':
                active_ij[0] = (active_ij[0] + 1) % 9
            elif s == 'left':
                active_ij[1] = (active_ij[1] + 9 - 1) % 9
 
            cursor.pos=vec(active_ij[1],active_ij[0],0)
            if len(nums[active_ij[0]][active_ij[1]].text) > 0:
                cursor.pos.x += font_w2
        elif s == 'esc':
            active_ij = [-1,-1]
            cursor.visible=False
        elif len(s)==1:
            #check valid entry only after clicking away
            s = int(s)
            if s > 0 and s <= 9:
                if len(the_num.text) == 0:
                    the_num.text += s # append new character
                    board[active_ij[0]][active_ij[1]] = int(s)
                    cursor.pos.x += font_w2
                    is_val = is_valid_at(active_ij[0],active_ij[1])
                    if not is_val:
                        the_num.color=invalid_color
                        invalid_warning.visible=True
                    else:
                        the_num.color=valid_input_num_color
                        invalid_warning.visible=False
        
scn.bind('keydown', key_input)

#need to make these input values fixed, check if visible
#only after solved do I make other numbers visible. If I solve with logic
#first, then make those visible before brute force


#eventually: try some logic to make it go faster than brute force

def increment_entry_at(row,col):
    #nums.text='' are non-fixed entries, increment only non-fixed entries
    global board
    if len(nums[row][col].text) == 0:
        if board[row][col] == 9:
            board[row][col] = init_entry_at(row,col)
        elif board[row][col] <= 0:
            board[row][col] = 1
        else:
            board[row][col] += 1
        return True
    else:
        return False

def is_filled():
    for i in range(9):
        for j in range(9):
            if board[i][j] <= 0:
                return False
                
    return True

def is_solved():
    return is_filled() and is_valid_board()

#solve game
def solve_board():
    #algorithm to solve if len(nums[i][j].text) == 0
    #goes across each row first: varies column more rapidly
    global nums, invalid_warning
    max_iters = 81000000   #why not 9**3? not sure, maybe it includes going back
    
    step = 1
    idx = 1 #from 1 to 82
    n = 0
    while idx >= 1 and n <= max_iters:
        n += 1
        row = int((idx - 1)/9)
        col = (idx - 1)%9
        is_incr = increment_entry_at(row,col)
        if is_incr:
            if board[row][col] > 0:
                if is_valid_at(row,col):
                    step = 1
                    idx += step
            else:
                step = -1
                idx += step
        else:
            idx += step
        
        if idx >= 82:
            if is_solved():
                for i in range(9):
                    for j in range(9):
                        if len(nums[i][j].text) == 0:
                            nums[i][j].color = num_color
                        nums[i][j].text = board[i][j]
            else:
                invalid_warning.visible=True
    
    if n >= max_iters:
        print('Hit max # of solutions:',max_iters)
                        
#button to run at bottom, make cursor invisible and active_ij=[-1,-1]
def solver_button(the_button):
    global ok_to_input, cursor
    cursor.visible=False
    if 'Solve' in the_button.text:
        the_button.text = 'Reload Page To Reset'
        ok_to_input = False
        solve_board() 
    
button(text=' Solve Sudoku Game ',pos=scn.caption_anchor,bind=solver_button)

