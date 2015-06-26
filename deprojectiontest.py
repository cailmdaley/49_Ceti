import numpy as np
from astropy.io import fits
from subprocess import call
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from pylab import *
import math
plt.close()

outer=100.0
inner=outer/2.5
pa=np.radians(-71.5)
inc=np.radians(79.1)
ratio=np.cos(inc)
major=100.0
minor=major*ratio
flux=6.0
dd=5.0 #disk diameter of 49 Ceti in degrees

imsize=301
gauss=np.zeros([imsize, imsize])
d=np.zeros([imsize, imsize])
xs=np.arange(imsize)-imsize/2.0
ys=np.arange(imsize)-imsize/2.0

#Disk
for i in range(0,imsize):
    for j in range(0,imsize):
        gauss[i,j]=flux*np.exp(-1.0*((xs[i]*np.cos(pa)-ys[j]*np.sin(pa))**2/(2.0*major**2) + (xs[i]*np.sin(pa) + ys[j]*np.cos(pa))**2/(2.0*minor**2)))
        d[i,j]= np.sqrt(((xs[i]*cos(pa) - ys[j]*np.sin(pa)))**2+((ys[j]*cos(pa) + xs[i]*sin(pa))/(minor/major)) **2)


#inner and outer radius
gauss[np.where((d<=inner) | (d>=outer))]=0
e=gauss/np.sum(gauss)*0.017

plt.imshow(e, origin='lower')
plt.savefig('virtual.pdf')
plt.show()


#
#Circular averaging
#drange = np.linspace(0, np.max(d))
#flux = np.zeros((len(drange)))
#for i  in range(len(drange)-1):
#     flux[i]=np.mean(e[np.where((d>=drange[i]) & (d<drange[i+1]))])
#plt.plot(drange,flux, '.g')
#plt.savefig('imcircavg.pdf')
#plt.show()



hdu=fits.PrimaryHDU(e)
hdu.header['CRVAL1']=15.0*(1.0+34.0/60.0+37.870/3600)
hdu.header['CRVAL2']=-15.0-40.0/60.0-34.94/3600.0
hdu.header['CDELT1']=dd/(-2.0*outer*3600.0)
hdu.header['CDELT2']=dd/(2.0*outer*3600.0)
hdu.header['CRPIX1']=(imsize+1)/2.0
hdu.header['CRPIX2']=(imsize+1)/2.0
hdu.header['EPOCH']=2000.0
hdu.writeto('49cetimodel.fits')
call("fits in=49cetimodel.fits op=xyin out=49cetimodel.im", shell=True)
call("uvmodel model=49cetimodel.im vis=49ceti.vis options=replace out=49cetimodel.vis", shell=True)
call("fits in=49cetimodel.vis op=uvout out=49cetimodel.uvf", shell=True)

