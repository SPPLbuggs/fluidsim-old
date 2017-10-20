import numpy as np
from scipy.io import FortranFile
import glob
import sys

x = np.fromfile('output/meshx.dat',dtype=float)
y = np.fromfile('output/meshy.dat',dtype=float)
t = np.fromfile('output/time.dat', dtype=float)

nx = len(x)
ny = max(len(y),1)
ts = len(t)

phi = np.zeros([nx,ny,ts])
ne = np.zeros([nx,ny,ts])
ni = np.zeros([nx,ny,ts])
nt = np.zeros([nx,ny,ts])
nm = np.zeros([nx,ny,ts])

temp = np.fromfile('output/phi.dat',dtype=float)
Ex = temp.reshape([ts, ny, nx])

temp = np.fromfile('output/ni.dat',dtype=float)
ni = temp.reshape([ts, ny, nx])

temp = np.fromfile('output/ne.dat',dtype=float)
ne = temp.reshape([ts, ny, nx])

temp = np.fromfile('output/nt.dat',dtype=float)
nt = temp.reshape([ts, ny, nx])
nt = nt / ne * 1.5

temp = np.fromfile('output/nm.dat',dtype=float)
nm = temp.reshape([ts, ny, nx])

np.savez('data/data.npz', x=x, y=y, t=t, phi=phi, ne=ne, ni=ni, nt=nt, nm=nm)
