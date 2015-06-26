import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import math
import pyfits
plt.close()

#diskimage = fits.open('60deg.uvf')
#header=pyfits.getheader('60deg.uvf')
diskimage = fits.open('49ceti.uvf')
header=pyfits.getheader('49ceti.uvf')
data = diskimage[0].data
rot=math.radians(30)
inc=math.acos(0.8)
freq = header['CRVAL4']
u = data.field(0)
u1 = u*freq*1e-03
v = data.field(1)
v1 = v*freq*1e-03
ww = data.field(2)
uu=math.cos(inc)*(u*math.cos(-rot)-v*math.sin(-rot))
vv=(u*math.sin(-rot)+v*math.cos(-rot))
rcw=data.field(5)
rlim=np.array(rcw)
rl=rlim[:,0,0,0,0,0]
im=rlim[:,0,0,0,0,1]
wt=rlim[:,0,0,0,0,2]
amp=(rl**2+im**2)**(0.5)
ubins=np.linspace(np.min(u1), np.max(u1), 20)
vbins=np.linspace(np.min(v1), np.max(v1), 20)
rlavg=np.zeros((len(ubins),len(vbins)))
imavg=np.zeros((len(ubins),len(vbins)))
p=np.zeros((len(ubins), len(vbins)))

#Figuring out WHAT IS GOING ON
#u1 = u[50000:50010]
#v1 = v[50000:50010]
#rl1 = rl[50000:50010]
#im1 = im[50000:50010]
#amp1 = amp[50000:50010]
#for i in range(len(ubins)-1):
#    for j in range(len(vbins)-1):
#        if np.size(np.where( (u1>=ubins[i]) & (u1<=ubins[i+1]) & (v1>=vbins[j]) & (v1<=vbins[j+1]) ))>0:
#            rlavg[i,j]=np.mean(rl1[np.where( (u1>=ubins[i]) & (u1<=ubins[i+1]) & (v1>=vbins[j]) & (v1<=vbins[j+1]) )] )
#            imavg[i,j]=np.mean(im1[np.where( (u1>=ubins[i]) & (u1<=ubins[i+1]) & (v1>=vbins[j]) & (v1<=vbins[j+1]) )] ) 
#        else:
#            rlavg[i,j]=0
#            imavg[i,j]=0
#ampavg=(rlavg**2+imavg**2)**(0.5)
##plt.pcolor(ubins,vbins, np.transpose(ampavg))
#plt.imshow(np.transpose(ampavg), origin='lower',  extent=[np.min(u), np.max(u), np.min(v), np.max(v)])
#plt.show()


#Scatter Plot
#plt.scatter(u,v,c=amp)
#plt.title('Amplitude as a function of (u,v)')
#plt.xlabel('u')
#plt.ylabel('v')
#v=[-2.0e-6,2.0e-6,-2.0e-06,2.0e-06]
#plt.axis(v)
#plt.savefig('amplitude_scatter_plot_60deg')
#plt.show()

#Array

for i in range(len(ubins)-1):
    for j in range(len(vbins)-1):
        if np.size(np.where( (u1>=ubins[i]) & (u1<=ubins[i+1]) & (v1>=vbins[j]) & (v1<=vbins[j+1]) ))>0:
            rlavg[i,j]=np.mean(rl[np.where( (u1>=ubins[i]) & (u1<=ubins[i+1]) & (v1>=vbins[j]) & (v1<=vbins[j+1]) )] )
            imavg[i,j]=np.mean(im[np.where( (u1>=ubins[i]) & (u1<=ubins[i+1]) & (v1>=vbins[j]) & (v1<=vbins[j+1]) )] ) 
        else:
            rlavg[i,j]=0
            imavg[i,j]=0
ampavg=(rlavg**2+imavg**2)**(0.5)

plt.imshow(np.transpose(ampavg), origin='lower',  extent=[np.min(u1), np.max(u1), np.min(v1), np.max(v1)])
plt.xlabel('u ($k\lambda$)')
plt.ylabel('v ($k\lambda$)') 
plt.savefig('49ceti visibilities')
plt.show()

#Circular Averaging
#d=(u**2+v**2)**0.5
#drange = np.linspace(0, np.max(d) )
#ampavg=np.zeros((len(drange)))
#for i in range(len(drange)-1):
#    if np.size(np.where( (d>=drange[i]) & (d<=drange[i+1]) ))>0:
#        ampavg[i]=np.mean(amp[np.where( (d>=drange[i]) & (d<=drange[i+1]) )])
#    else: ampavg[i]=0
#plt.plot(drange,ampavg, '.g')
#plt.savefig('circularaveraging.pdf')
#plt.show()

