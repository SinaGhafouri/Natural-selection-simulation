import random as rn
import numpy as np
import matplotlib.pyplot as plt
import time

a = 50  #Size of the land = n*n (Not considering the edges).
m = 50  #Indicates amount of resources.
day = 10  #Number of days.
E = 100  #Energy of the sim.
N = 50   #Number of sims.
v1, v1_percentage, v2 = 1, 50, 2    #Velocity quntities.
s1, s1_percentage, s2 = 0, 50, 3   #Sensability quantities.

def vel(v1, v1_percentage, v2):     #Velocity of the sims.
    V1 = int(np.floor(N/100*v1_percentage))*[int(v1)]
    V2 = int(np.ceil(N/100*(100-v1_percentage)))*[int(v2)]
    return V1+V2

def Sense_Quantities(s1, s1_percentage, s2):     #Sensability of the sims.
    S1 = int(np.floor(N/100*s1_percentage))*[int(s1)]
    S2 = int(np.ceil(N/100*(100-s1_percentage)))*[int(s2)]
    return S1+S2

def sense(S, x, y, ran):  #Senseability function in random walk.
    for s in range(S, 0, -1):
        if L[d][int((x+s)%(a+1))][y] == 1:
            ran = 0
        elif L[d][x][int((y+s)%(a+1))] == 1:
            ran = 1
        elif L[d][int((x-s)%(a+1))][y] == 1:
            ran = 2
        elif L[d][x][int((y-s)%(a+1))] == 1:
            ran = 3
        elif L[d][int((x+s)%(a+1))][int((y+s)%(a+1))] == 1:
            ran = 4
        elif L[d][int((x-s)%(a+1))][int((y-s)%(a+1))] == 1:
            ran = 5
        elif L[d][int((x-s)%(a+1))][int((y+s)%(a+1))] == 1:
            ran = 6
        elif L[d][int((x+s)%(a+1))][int((y-s)%(a+1))] == 1:
            ran = 7
        else: ran = ran
    return ran

def sim(N):  #This function creats a sim and give it a coordinate on the edge.
    global rSim
    rSim = [[], []]    #Cordinate of sims.
    for _ in range(N):
        ran = np.floor(rn.random()*4)
        if ran == 0:
            rSim[0].append(0)
            rSim[1].append(int(rn.random()*(a+2)))
        elif ran == 1:
            rSim[0].append(int(rn.random()*(a+2)))
            rSim[1].append(0)
        elif ran == 2:
            rSim[0].append(a+2-1)
            rSim[1].append(int(rn.random()*(a+2)))
        elif ran == 3:
            rSim[0].append(int(rn.random()*(a+2)))
            rSim[1].append(a+2-1)
    return rSim     #rSim[x or y][sim's tag]

def land(day, a, m):
    global L, xRe, yRe
    L = np.zeros((day, a+2, a+2))
    for j in range(day):
        xRe = []
        yRe = []
        i = 0
        while i < m:
            new_x = int(np.floor(rn.random()*(a)))
            new_y = int(np.floor(rn.random()*(a)))
            for k in range(len(xRe)):
                if new_x == xRe[k] and new_y == yRe[k]:
                    break
            else:
                L[j][new_x+1][new_y+1] = 1
                xRe.append(new_x)
                yRe.append(new_y)
                i += 1
    return L
 
t1 = time.time()
x0n0, y0n0 = sim(N)[0], sim(N)[1] #First day's posiiton of each sim.
x0n, y0n = [], []

for n in range(N):
    x0n.append(x0n0[n])
    y0n.append(y0n0[n])
        
land(day, a, m)
path_test = -1    #For not going back (It can be any number except 0,1 and 2).
EC = np.zeros((day, N))     #Energy consumption.
dead_count_by_day = np.zeros(day)
dead_count_by_vel_and_day = np.zeros((day, 2))    #2 different velocities.
dead_count_by_sense_and_day = np.zeros((day, 2))    #2 different sensabilities.
dead_id = []    #Those sims which are dead.
edge_id = []    #Those sims which are on the edge.
Xn = []
Yn = []
V = vel(v1, v1_percentage, v2)  #We will use this in the loop!
V_test = vel(v1, v1_percentage, v2) #For dead count.
S = Sense_Quantities(s1, s1_percentage, s2)  #We will use this in the loop!
S_test = Sense_Quantities(s1, s1_percentage, s2) #For dead count.
for d in range(day):
    #print('\nNew day:', d+1)
    dead_id.append([])
    edge_id.append([])
    Xn.append([])
    Yn.append([])
    if d != 0:
        for nn in range(N):
            Xn[d].append([Xn[d-1][nn][-1]])
            Yn[d].append([Yn[d-1][nn][-1]])
    else:   #For the first day.
        [Xn[d].append([x0n[i]]) for i in range(N)]
        [Yn[d].append([y0n[i]]) for i in range(N)]
    PT = N*[-1] #(Path Test) For not going back (It can be any number except 0,1 and 2(random numbers for out latter random walk)).
    for e in range(E):
        v = np.array(V)
        v_test = np.array(V_test)
        s = np.array(S)
        s_test = np.array(S_test)
        n = 0
        h=0
        while True: #One step random walk for the nth sim.
            if n >= N: break
            if n in dead_id[d]:
                n += 1
                continue
            if n in edge_id[d]:
                n += 1
                continue

            r = int(rn.random()*8)
            r = sense(s[n], Xn[d][n][e], Yn[d][n][e], r)
            
            if r == 0:
                #if r == int((PT[n]+2)%8): continue
                Xn[d][n].append((Xn[d][n][e] + 1)%(a+1))
                Yn[d][n].append(Yn[d][n][e])
                PT[n] = r

            elif r == 1:
                #if r == int((PT[n]+2)%8): continue
                Xn[d][n].append(Xn[d][n][e])
                Yn[d][n].append((Yn[d][n][e] + 1)%(a+1))
                PT[n] = r

            elif r == 2:
                #if r == int((PT[n]+2)%8): continue
                Xn[d][n].append((Xn[d][n][e] - 1)%(a+1))
                Yn[d][n].append(Yn[d][n][e])
                PT[n] = r

            elif r == 3:
                #if r ==  int((PT[n]+2)%8): continue
                Xn[d][n].append(Xn[d][n][e])
                Yn[d][n].append((Yn[d][n][e] - 1)%(a+1))
                PT[n] = r            
            elif r == 4:
                #if r == int((PT[n]+2)%8): continue
                Xn[d][n].append((Xn[d][n][e] + 1)%(a+1))
                Yn[d][n].append((Yn[d][n][e] + 1)%(a+1))
                PT[n] = r

            elif r == 5:
                #if r == int((PT[n]+2)%8): continue
                Xn[d][n].append((Xn[d][n][e] - 1)%(a+1))
                Yn[d][n].append((Yn[d][n][e] - 1)%(a+1))
                PT[n] = r

            elif r == 6:
                #if r == int((PT[n]+2)%8): continue
                Xn[d][n].append((Xn[d][n][e] - 1)%(a+1))
                Yn[d][n].append((Yn[d][n][e] + 1)%(a+1))
                PT[n] = r

            elif r == 7:
                #if r ==  int((PT[n]+2)%8): continue
                Xn[d][n].append((Xn[d][n][e] + 1)%(a+1))
                Yn[d][n].append((Yn[d][n][e] - 1)%(a+1))
                PT[n] = r

            if L[d][Xn[d][n][e+1]][Yn[d][n][e+1]] == 1:   #this sim survived today.
                #print('*** Day',d+1,' sim',n+1, 'Hoorays ^_^ ***')
                edge_id[d].append(n)
                L[d][Xn[d][n][e+1]][Yn[d][n][e+1]] = 0
                ###To send them back to the nearest edge:###
                dis = [abs(1-Xn[d][n][e+1]),abs(a+1-Xn[d][n][e+1]),abs(1-Yn[d][n][e+1]),abs(a+1-Yn[d][n][e+1])]
                Min_dis = dis.index(min(dis))
                if Min_dis == 0:
                    Xn[d][n].append(0)
                    Yn[d][n].append(Yn[d][n][e+1])
                if Min_dis == 1:
                    Xn[d][n].append(a+1)
                    Yn[d][n].append(Yn[d][n][e+1])
                if Min_dis == 2:
                    Xn[d][n].append(Xn[d][n][e+1])
                    Yn[d][n].append(0)
                if Min_dis == 3:
                    Xn[d][n].append(Xn[d][n][e+1])
                    Yn[d][n].append(a+1)

            else: EC[d][n] += 1*v[n]**2 + s[n]    #K is proportional to v^2. Sensability take energy too!
            if EC[d][n] >= E:
                ###To send them back to the nearest edge:###
                dis = [abs(1-Xn[d][n][e+1]),abs(a+1-Xn[d][n][e+1]),abs(1-Yn[d][n][e+1]),abs(a+1-Yn[d][n][e+1])]
                Min_dis = dis.index(min(dis))
                if Min_dis == 0:
                    Xn[d][n].append(0)
                    Yn[d][n].append(Yn[d][n][e+1])
                if Min_dis == 1:
                    Xn[d][n].append(a+1)
                    Yn[d][n].append(Yn[d][n][e+1])
                if Min_dis == 2:
                    Xn[d][n].append(Xn[d][n][e+1])
                    Yn[d][n].append(0)
                if Min_dis == 3:
                    Xn[d][n].append(Xn[d][n][e+1])
                    Yn[d][n].append(a+1)
                #print('*** Day',d+1,' sim',n+1, 'is Dead X_X ***')
                dead_count_by_day[d] += 1
                dead_id[d].append(n)
                if n < int(v1_percentage*N/100): dead_count_by_vel_and_day[d][0]+=1
                elif v_test[n] == v2: dead_count_by_vel_and_day[d][1]+=1
                if n < int(s1_percentage*N/100): dead_count_by_sense_and_day[d][0]+=1
                else: dead_count_by_sense_and_day[d][1]+=1
            
            v[n] -= 1
            if v[n] > 0:
                continue
            
            n += 1  #The while loop counter.

t2 = time.time()
print('\nTime elapsed = ' ,"%.2f" % (t2-t1), 'Sec')

print('\nDead count by day:\n', dead_count_by_day)
print('\nDead count by velocity and day [v1 =',v1,', v2 =',v2,'] and v1/v2 ratio = ',v1_percentage,'/',(100-v1_percentage),':\n', dead_count_by_vel_and_day)
print('\nDead count by sensability and day [s1 =',s1,', s2 =',s2,'] and s1/s2 ratio = ',s1_percentage,'/',(100-s1_percentage),':\n', dead_count_by_sense_and_day)


import matplotlib.animation as animation
fig = plt.figure()
ax = plt.axes(xlim=(0, a+2), ylim=(0, a+2))
line_land, = ax.plot([], [], 'g.', linestyle="")
lines = []
for i in range(N):
    l = ax.plot([], [], 'o')[0]
    lines.append(l)

simx = []
simy = []
for i in range(len(Xn[0])):
    simx.append([])
    simy.append([])
    for j in range(len(Xn[0][i])):
        simx[i].append(Xn[0][i][j])
        simy[i].append(Yn[0][i][j])
Lx = []
Ly = []
for i in range(a+1):
    for j in range(a+1):
        if L[0][i][j] == 1:
            Lx.append(i)
            Ly.append(j)
            
X, Y = [], []
def animate(i):
    for j in range(N):
        try:
            X.append([])
            Y.append([])
            
            xs = simx[j][i]
            ys = simy[j][i]
            
            X[j].append(xs)
            Y[j].append(ys)
        except: continue    
    xlist = [X[j][-1] for j in range(N)]
    ylist = [Y[j][-1] for j in range(N)]

    for lnum,line in enumerate(lines):
        line.set_data(xlist[lnum], ylist[lnum])
        
    LANDx = Lx
    LANDy = Ly
    
    line_land.set_data(LANDx, LANDy)
    return line, line_land,
ani = animation.FuncAnimation(fig, animate, interval=100)
plt.show()
