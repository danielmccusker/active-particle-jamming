import numpy as np
import pylab as py
import matplotlib.pyplot as plt
import matplotlib.collections
import csv
import os
from matplotlib import animation
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from PIL import Image

np.set_printoptions(threshold=np.nan)

num = 512
path = "positions2.dat"
infile = open(path, "rb")
reader = csv.reader(infile, delimiter = "\t",quoting=csv.QUOTE_NONNUMERIC)
x = []
y = []
r = []
for row in reader:
    x.append(row[0])
    y.append(row[1])
    r.append(row[2])
infile.close()

fig, ax = plt.subplots()

xtemp = x[:num]
ytemp = y[:num]
rtemp = r[:num]
xy = np.column_stack((xtemp,ytemp))
patches = []
patches = [plt.Circle(center, size) for center, size in zip(xy, rtemp)]
coll = matplotlib.collections.PatchCollection(patches, animated=False, facecolors='black')
ax.add_collection(coll)
ax.autoscale_view(True)
ax.set_aspect('equal')
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
fig.set_size_inches([5,5])
py.tight_layout()

def init():
    return coll,

def animate(i):
    patches = []
    xtemp = []
    ytemp = []
    rtemp= []
    xtemp = x[num*i:num*(i+1)]
    ytemp = y[num*i:num*(i+1)]
    rtemp = r[num*i:num*(i+1)]
    xy = np.column_stack((xtemp,ytemp))
    patches = [plt.Circle(center, size) for center, size in zip(xy, rtemp)]
    coll = matplotlib.collections.PatchCollection(patches, animated=True, facecolors='black')
    ax.add_collection(coll)
    
    return coll,

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=len(x)/num, interval=1e1, blit=True, repeat=False)
anim.save("videoTest.mp4")

fig.canvas.draw()
w,h = fig.canvas.get_width_height()
buf = np.fromstring (fig.canvas.tostring_argb(), dtype=np.uint8)
buf.shape = (w,h,4)
#canvas.tostring_argb give pixmap in ARGB mode. Roll the ALPHA channel to have it in RGBA mode
buf = np.roll (buf,3,axis = 2)
#print buf

ft = np.fft.fft2(buf)
ft = np.fft.fftshift(abs(ft))
#print ft
#plt.plot(ft)

w, h, d = buf.shape
im = Image.frombytes( "RGBA", ( w ,h ), buf.tostring( ) )
img = np.fft.fftshift(abs(np.fft.fft2(im)))
py.clf()
py.imshow(img)
py.show()
#plt.show()
