import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
geom = []
temp = []

fp = open("geom","r+")
for i, line in enumerate(fp):
    if i == 0:
        geom.append(line.split())
fp.close()

#getting data
rows = geom[0][0]
writeInterval = int(geom[0][1])
Length = float(geom[0][2])
Width = float(geom[0][3])
N = int(geom[0][4])
M = int(geom[0][5])
times = int(geom[0][6])

x = np.linspace(0, Length, N)
y = np.linspace(0, Width, M)
[X, Y] = np.meshgrid(x,y)

data = np.genfromtxt("T", delimiter=",")
minTemp = data.min()
maxTemp = data.max()
if (minTemp==maxTemp):
    maxTemp = maxTemp+ 1/maxTemp


def getTempForTime(i):
    return data[M*i:M*(i+1),:]
    #return np.flip(data[M*i:M*(i+1),:],axis=0)

figHeight = int(Width*10/(Width+Length))
figLength = int(Length*10/(Width+Length))

fig = plt.figure(figsize=(figLength,figHeight))
ax = plt.axes(xlim=(0, Length), ylim=(0, Width), xlabel='Length', ylabel='Width')

cvals = np.linspace(minTemp,maxTemp,20)      # set contour values 
cont = plt.contourf(x, y, getTempForTime(0), cvals)    # first image on screen
plt.colorbar()
plt.savefig(str(times*writeInterval).zfill(5) )

# animation function
def animate(i):
    global cont
    T = getTempForTime(i)
    for c in cont.collections:
        c.remove()  # removes only the contours, leaves the rest intact
    cont = plt.contourf(x, y, T, cvals)
    plt.title('t = %i:  %.2f' % (i*writeInterval,times*writeInterval))
    return cont

anim = animation.FuncAnimation(fig, animate, frames=times, repeat=False)
anim.save('animation.mp4', writer=animation.FFMpegWriter())

