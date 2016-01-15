import numpy as np
import matplotlib.pyplot as plt
from pylab import *
import math
plt.clf()


e=np.zeros([251,251])
d=np.zeros([251,251])
#p=np.zeros([110,110])
#d2=np.zeros([110,110])
e1=np.arange(0, 100)
d1=np.arange(0, 100)
def degToRad(i):
    return (math.pi/180)*i

for x in range(0,251):
    for y in range(0,251):
        if ((x-125.0)*np.cos(120)+(y-125)*np.sin(120))**2/100.0+((x-125.0)*np.sin(120)-(y-125)*np.cos(120))**2/25.0>=1 and (((x-125.0)*np.cos(120)+(y-125)*np.sin(120))**2/10000.0+((x-125.0)*np.sin(120)-(y-125)*np.cos(120))**2/2500.0)<=1:
        #d[x,y]=np.sqrt(((1/np.cos(degToRad(60)))*(y-125))**2+(x-125)**2)
        #if d[x,y]>=10 and d[x,y]<=100:
        #if (x-125.0)**2/100.0+(y-125)**2/25.0>=1 and (x-125.0)**2/10000.0+(y-125)**2/2500.0<=1:
            #e[x,y]=2.0*np.exp(-1.0*((x-125.0)**2/(2.0*100.0**2)+(y-125.0)**2/(2.0*50.0**2)))
            e[x,y]=2.0*np.exp(-1.0*(((x-125.0)*np.cos(120)+(y-125)*np.sin(120))**2/(2.0*100.0**2)+((x-125.0)*np.sin(120)-(y-125)*np.cos(120))**2/(2.0*50.0**2)))
        else: e[x,y]=0

    
#e1[d1]=2.0*np.exp(-1.0*(d1-125.0)**2/(2.0*100.0**2))
#plt.plot(d1,e1)
#for i in range(110):
    #p[i]=np.mean(e[np.where((d<i) & (d>=(i-1)))])
    #d2[i]=i-0.5

#print np.sum(e[np.where((d<9))])




#plt.plot(d2,p,'.m')
plt.title('Disk Image')
plt.xlabel('AU')
plt.ylabel('AU')
plt.imshow(e)
plt.savefig('diskimage.pdf')
plt.show()

