import numpy as np
import matplotlib.pyplot as plt
from matplotlib.image import NonUniformImage
from matplotlib import cm
from scipy.io import FortranFile
import glob
from matplotlib import animation
#from matplotlib import rcParams
#rcParams.update({'figure.autolayout': True})

size = 12
med_size = 13
big_size = 14

plt.rc('font', size=size)
plt.rc('axes', titlesize=size)
plt.rc('axes', labelsize=med_size)
plt.rc('xtick', labelsize=size)
plt.rc('ytick', labelsize=size)
plt.rc('legend', fontsize=size)
plt.rc('figure', titlesize=big_size)
plt.rcParams['figure.figsize'] = (4.5, 3)
#plt.rcParams['figure.autolayout'] = True

cm_subsection = np.linspace(0.0, 1.0, 4) 
colors = [ cm.viridis(x) for x in cm_subsection ]

path = 'output/'
x = np.fromfile('output/meshx.dat',dtype=float)
y = np.fromfile('output/meshy.dat',dtype=float)
t = np.fromfile('output/time.dat', dtype=float)

nx = len(x)
ny = len(y)
ts = len(t)

Ex = np.zeros([nx,ny,ts])
Ey = np.zeros([nx,ny,ts])
ne = np.zeros([nx,ny,ts])
ni = np.zeros([nx,ny,ts])
nt = np.zeros([nx,ny,ts])
nm = np.zeros([nx,ny,ts])

temp = np.fromfile('output/ne.dat',dtype=float)
Ex = temp.reshape([ts, ny, nx])

interp = 'bilinear'
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,4))
im = NonUniformImage(ax, interpolation=interp, cmap='plasma', 
                     extent=(x[0], x[-1], y[0], y[-1]))
im.set_data(x, y, Ex[0,:,:])
ax.images.append(im)
ax.set_xlim(x[0]-0.1, x[-1]+0.1)
ax.set_ylim(y[0]-0.1, y[-1]+0.1)
tx = plt.title(r'$n_e$ : time = 0 $\mu$s')
ax.set_xlabel('x [mm]')
ax.set_ylabel('y [mm]')
plt.tight_layout()

def anim(i):
    im.set_data(x, y, Ex[i,:,:])
    ax.set_xlim(x[0]-0.1, x[-1]+0.1)
    ax.set_ylim(y[0]-0.1, y[-1]+0.1)
    tx.set_text(r'$n_e$ : time = {:.2f} $\mu$s'.format(t[i]))
    im.autoscale()

ani = animation.FuncAnimation(fig, anim, frames = ts, interval = 75)

plt.show()
