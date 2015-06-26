import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import math
import pyfits
from subprocess import call
import os
plt.close()


#INPUT VARIABLES
#print "This code will read in a FITS file and deproject its visibilities. It can handle line emission and continuum data, and can radially average the visibilities to create a plot of flux as a function of radially deprojected distance and/or image the disk with a single image or multiple channel maps."  
#fileinp = raw_input('Please enter the filename: ')
#print
#print "If there is an offset between the center of the source and the center of the image, please enter the shift required to center the source. (if the source center is at (-1.0,-1.0) you would enter (1.0,1.0)"
#print
#deltara = math.radians( input('RA shift in arcsecs: ')/3600.0)
#deltadec = math.radians( input('Dec shift in arcsecs: ')/3600.0) 
#pa = math.radians( input('Please enter the SKY PA in degrees: '))
#inc = math.radians( input('Please enter the inclination of the disk in degrees: '))
#datatype = raw_input('Please enter \'line\' for line emission data or \'cont\' for continuum data: ')
#outputtype = raw_input('Please enter \'rad\' for radial averaging, \'im\' for imaging, or \'both\' for both: ')
#if outputtype in ('im', 'both'):
#    image=raw_input('Please enter the filename you would like your image(s) to have, including the .uvf extension: ')
#if outputtype in ('rad', 'both'):
#    bins = input('Please enter the number of bins you would like: ')
#    binning = raw_input('Please enter \'log\' for logarithmic binning or \'lin\' for linear binning: ')
#    try:
#        xlim = input('If you would like to set a value for the maximum distance plotted, please enter it in $k\lambda$ here :')
#    except SyntaxError:
#        xlim = None
#model = raw_input('If you would like to compare with a model, please include its filename here: ')
#if len(model)>4:
#   modpa = math.radians( input('Please enter the SKY PA of the model in degrees (positive angles only): '))
#   modinc = math.radians( input('Please enter the inclination of the model in degrees: ')
fileinp='49ceti.uvf'
outputtype='rad'
bins=20
binning='lin'
spa=108.5
pa=math.radians(spa)
inc=math.radians(79.1)
model='.uvf'
modpa=math.radians(108.5)
modinc=math.radians(79.1)
image='age'
try:
    xlim = input('If you would like to set a value for the maximum distance plotted, please enter it in $k\lambda$ here :')
except SyntaxError:
    xlim = None

deltara=0
deltadec=0

#Reading in fits file/data (assuming ALMA data) and converting to kilolambda
dfits = fits.open(fileinp)
data = dfits[0].data
vis = dfits[0].data['data']
header = dfits[0].header
freq = header['CRVAL4']
u = data['UU']
v = data['VV']
#chans = vis[0,0,0,
rlxx = (np.ravel(vis[:,0,0,0,0,0], order='F'))*1e03 
rlyy = (np.ravel(vis[:,0,0,0,1,0], order='F'))*1e03 
imxx = (np.ravel(vis[:,0,0,0,0,1], order='F'))*1e03
imyy = (np.ravel(vis[:,0,0,0,1,1], order='F'))*1e03
rl = (rlxx + rlyy)/2.0
im = (imxx + imyy)/2.0
wt = np.ravel(vis[:,0,0,0,0,2], order='F')
amp=np.sqrt(rl**2+im**2)

#Calculating phase shifts due to position offset
deltaphase = 2*math.pi*(u*deltara + v*deltadec)
phase = np.arctan(im/rl)+deltaphase
rlcor = amp*np.cos(phase)
imcor = amp*np.sin(phase)

#Rotation of (u,v) coordinates back to PA=0 and deprojection
u1 = (u*np.cos(pa)-v*np.sin(pa))*np.cos(inc)
v1 = (u*np.sin(pa)+v*np.cos(pa))
#Rotation of (u,v) coordinates back to sky PA 
u2 = u1*np.cos(-1*pa)-v1*np.sin(-1*pa)
v2 = u1*np.sin(-1*pa)+v1*np.cos(-1*pa)
uu = u2*freq*1e-03
vv = v2*freq*1e-03
depdist = (uu**2 + vv**2)**0.5

#Radial Averaging
if outputtype in ('rad', 'both'):
    if xlim:
        if binning == 'lin':
            distrange = np.linspace(np.min(depdist), xlim, bins)
        if binning == 'log':
            distrange = np.logspace(np.log10(np.min(depdist)), np.log10(xlim), bins)
    else:
        if binning == 'lin':
            distrange = np.linspace(np.min(depdist), np.max(depdist), bins)
        if binning == 'log':
            distrange = np.logspace(np.log10(np.min(depdist)), np.log10(np.max(depdist)), bins)
    rlavg = np.zeros(bins)
    imavg = np.zeros(bins)
    N = np.zeros(bins)
    for i in range(bins-1):
        if np.size(np.where( (depdist>=distrange[i]) & (depdist<=distrange[i+1]) ))==0:
            rlavg[i]=0
            imavg[i]=0
        else:
            rlavg[i] = np.sum((rlcor*wt)[np.where( (depdist>=distrange[i]) & (depdist<=distrange[i+1]) )]) / (np.sum(wt[np.where( (depdist>=distrange[i]) & (depdist<=distrange[i+1]) )]))

            imavg[i] = np.sum((imcor*wt)[np.where( (depdist>=distrange[i]) & (depdist<=distrange[i+1]) )]) / (np.sum(wt[np.where( (depdist>=distrange[i]) & (depdist<=distrange[i+1]) )]))
    ampavg = np.sqrt(rlavg**2 + imavg**2)


    #Standard Error of the Mean
    rlsamplestandev = np.zeros(bins)
    imsamplestandev = np.zeros(bins)
    for i in range(bins):
        if i<np.max(range(bins)):  
            rlsamplestandev[i] = np.std(rlcor[np.where( (depdist>=distrange[i]) & (depdist<=distrange[i+1]) )], ddof=1)
            imsamplestandev[i] = np.std(imcor[np.where( (depdist>=distrange[i]) & (depdist<=distrange[i+1]) )], ddof=1)
            N[i]=np.size(np.where( (depdist>=distrange[i]) & (depdist<=distrange[i+1]) ))
        else:
            rlsamplestandev[i] = np.std(rlcor[np.where( (depdist>=distrange[i]) & (depdist<=np.max(distrange)) )], ddof=1)
            imsamplestandev[i] = np.std(imcor[np.where( (depdist>=distrange[i]) & (depdist<=np.max(distrange)) )], ddof=1)
            N[i]=np.size(np.where( (depdist>=distrange[i]) & (depdist<=np.max(distrange)) ))
    rlSEM=rlsamplestandev/np.sqrt(N)
    imSEM=imsamplestandev/np.sqrt(N)
    
#Model
if len(model)>4:
    modfits = fits.open(model)
    moddata = modfits[0].data
    modvis = modfits[0].data['data']
    modheader = modfits[0].header
    modu = moddata['UU']
    modv = moddata['VV']
    modrl = (np.ravel(modvis[:,0,0,0,:,0], order='F'))*1e03 
    modim = (np.ravel(modvis[:,0,0,0,:,1], order='F'))*1e03
    modwt = np.ravel(modvis[:,0,0,0,:,2], order='F')
    modfreq = modheader['CRVAL4']
    modamp = np.sqrt((modrl**2+modim**2))

    modu1 = (modu*np.cos(modpa)-modv*np.sin(modpa))*np.cos(modinc)
    modv1 = modu*np.sin(modpa)+modv*np.cos(modpa)
    modu2 = modu1*np.cos(-1*modpa)-modv1*np.sin(-1*modpa)
    modv2 = modu1*np.sin(-1*modpa)+modv1*np.cos(-1*modpa)
    moduu = modu2*modfreq*1e-03
    modvv = modv2*modfreq*1e-03
    moddepdist = (moduu**2 + modvv**2)**0.5


    if binning == 'lin':
        moddistrange = np.linspace(np.min(moddepdist), np.max(moddepdist), bins)
    if binning == 'log':
        moddistrange = np.logspace(np.log10(np.min(moddepdist)), np.log10(np.max(moddepdist)), bins)
    modrlavg = np.zeros(bins)
    modimavg = np.zeros(bins)
    for i in range(bins-1):
        if np.size(np.where( (moddepdist>=moddistrange[i]) & (moddepdist<=moddistrange[i+1]) ))==0:
            modrlavg[i]=0
            modimavg[i]=0
        else:
            modrlavg[i] = np.sum((modrl*modwt)[np.where( (moddepdist>=moddistrange[i]) & (moddepdist<=moddistrange[i+1]) )]) / (np.sum(modwt[np.where( (moddepdist>=moddistrange[i]) & (moddepdist<=moddistrange[i+1]) )]))

            modimavg[i] = np.sum((modim*modwt)[np.where( (moddepdist>=moddistrange[i]) & (moddepdist<=moddistrange[i+1]) )]) / (np.sum(modwt[np.where( (moddepdist>=moddistrange[i]) & (moddepdist<=moddistrange[i+1]) )]))
    modampavg = np.sqrt(modrlavg**2 + modimavg**2)


#Plotting
plt.subplot(211)
plt.errorbar( distrange, rlavg, yerr=rlSEM, marker='.', fmt=' ' )
if len(model)>4:
    plt.plot(moddistrange, modrlavg) 
plt.title('Flux as a Function of Deprojected Distance')
plt.xlabel('$k\lambda$')
plt.ylabel('mJy')

plt.subplot(212)
plt.errorbar( distrange, imavg, yerr=imSEM, marker='.', fmt=' ' )
if len(model)>4:
    plt.plot(moddistrange, modimavg) 
plt.xlabel('$k\lambda$')
plt.ylabel('mJy')
#plt.savefig('49cetibinned.pdf')
plt.show()

#Creating visibility fits file
if len(image)>=4:
    data['UU'] = u2
    data['VV'] = v2
    dfits.writeto('deprojected_disk_image.uvf')
    call("fits in=deprojected_disk_image.uvf op=uvin out=deprojected_disk_image.vis", shell=True)
    call("invert vis=deprojected_disk_image.vis map=deprojected_disk_image.mp beam=deprojected_disk_image.bm cell=0.3 imsize=256 options=systemp,mfs robust=2", shell=True)
    call("clean map=deprojected_disk_image.mp beam=deprojected_disk_image.bm out=deprojected_disk_image.cl niters=300", shell=True)
    call("restor map=deprojected_disk_image.mp beam=deprojected_disk_image.bm model=deprojected_disk_image.cl out=deprojected_disk_image.cm", shell=True)
    call("cgdisp in=deprojected_disk_image.cm device=/xs labtyp=arcsec beamtyp=b,l,4", shell=True)

