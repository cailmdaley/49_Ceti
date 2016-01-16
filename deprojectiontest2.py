import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from pylab import *
import math
plt.close()

outer=240.0
inner=60.0
major=200.0
minor=100.0
flux=6.0
pa=np.radians(0)
inc=np.radians(0)
dd=5.0 #disk diameter of 49 Ceti in degrees

imsize=601
e=np.zeros([imsize, imsize])
d=np.zeros([imsize, imsize])
xs=np.arange(imsize)-imsize/2.0
ys=np.arange(imsize)-imsize/2.0

#Disk
for i in range(0,imsize):
    for j in range(0,imsize):
        e[i,j]= flux*np.exp(-1.0*((xs[i]*np.cos(pa)+ys[j]*np.sin(pa))**2/(2.0*major**2) + (xs[i]*np.sin(pa) - ys[j]*np.cos(pa))**2/(2.0*major**2)))
        d[i,j]= np.sqrt(((xs[i]*cos(pa) +ys[j]*np.sin(pa)))**2+(ys[j]*cos(pa) - xs[i]*sin(pa))**2)

e=np.transpose(e)
d=np.transpose(d)


#inner and outer radius
e[np.where((d<=inner) | (d>=outer))]=0
plt.imshow(e)
plt.xticks([])
plt.yticks([])
plt.savefig("/home/cdaley/Desktop/circle.png")
plt.show()

