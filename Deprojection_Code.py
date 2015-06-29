import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import math
import pyfits
from subprocess import call
import os
plt.close()

#INPUT VARIABLES
#print "This code will read in a FITS file and deproject its rlimwtibilities. It can handle line emission and continuum data, and can radially average the visibilities to create a plot of flux as a function of radially deprojected distance and/or image the disk with a single image or multiple channel maps."  
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
#   image = raw_input('Please enter the filename you would like your image(s) to have, not including the filetype: ')
#   if datatype == 'line':
#      imchans = 'line=channel,' + raw_input('Please enter the channels you would like to image in standard line miriad format, of the form  \'number of channels,starting channel index,channel width,step size.\' Indexing starts at 1: ')    
#       chandims = 'nxy=' + raw_input('Please enter the dimensions of the channel maps to be shown on each page, in the form \'columns,rows\': ')
#if outputtype in ('rad', 'both'):
#   if dataype == 'line':
#       chansin =input('Please enter the indices of the channel(s) you would like to average, keeping in mind that indexing starts at 0 rather than one. For example, one might enter \'5,12\' if they wanted channels 6 through 12. If you would only like to average one channnel, please enter the index of the channel and then the index of the next channel (e.g. \'2,3\') ')
#       chans = range(chans[0],chans[1])
#    bins = input('Please enter the number of bins you would like: ')
#    binning = raw_input('Please enter \'log\' for logarithmic binning or \'lin\' for linear binning: ')
#    try:
#        xlim = input('If you would like to set a value for the maximum distance plotted, please enter it in $k\lambda$ here: ')
#    except SyntaxError:
#        xlim = None
#model = raw_input('If you would like to compare with a model, please include its filename here: ')
#if len(model)>4:
#   modpa = math.radians( input('Please enter the SKY PA of the model in degrees (positive angles only): '))
#   modinc = math.radians( input('Please enter the inclination of the model in degrees: ')
fileinp='ceti250ch.uvf'
outputtype='im'
bins=15
binning='log'
spa=108.5
pa=math.radians(spa)
inc=math.radians(79.1)
model='.uvf'
modpa=math.radians(108.5)
modinc=math.radians(79.1)
image='49cetichannelmap'
deltara=0
deltadec=0
xlim = 300
datatype= 'line'
imchans = 'line=channel,'+'100,1,1,1'
chandims = 'nxy=6,3'
    


#Reading in fits file/data (assuming ALMA data)
dfits = fits.open(fileinp)
data = dfits[0].data
rlimwt = dfits[0].data['data']
header = dfits[0].header
freq = header['CRVAL4']
u = data['UU']
v = data['VV']
visnum=len(u)


#Reading in line emission visibilitites//weights and calculating phase shifts
if outputtype in ('rad', 'both'):
    if datatype == 'line':
        chans = range(0,250)
        rlxx = rlimwt[:,0,0,0,chans,0,0]
        rlyy = rlimwt[:,0,0,0,chans,1,0]
        imxx = rlimwt[:,0,0,0,chans,0,1]
        imyy = rlimwt[:,0,0,0,chans,1,1]
        wtxx = rlimwt[:,0,0,0,chans,0,2]
        wtyy = rlimwt[:,0,0,0,chans,1,2]
        rl = np.array((rlxx*wtxx+rlyy*wtyy)/(wtxx+wtyy))
        im = np.array((imxx*wtxx+imyy*wtyy)/(wtxx+wtyy))
        wt = 4.0/(1.0/wtxx + 1.0/wtyy)
        amp=np.sqrt(rl**2+im**2)

        deltaphase = 2*math.pi*(u*deltara + v*deltadec)
        deltaphase = (np.tile(deltaphase, (len(chans), 1))).T
        phase = np.arctan(im/rl)+deltaphase
        rlcor = np.cos(phase)*amp
        rlcor[np.where( rl<0 )] = -1.0*rlcor[np.where( rl<0 )]
        imcor = np.sin(phase)*amp
        imcor[np.where( rl<0)] = -1.0*imcor[np.where( rl<0 )]


    #Reading in continuum visibilitites//weights and calculating phase shifts
    if datatype == 'cont':
        chans = header['NAXIS4']
        rlxx = rlimwt[:,0,0,0,0,0]
        rlyy = rlimwt[:,0,0,0,1,0]
        imxx = rlimwt[:,0,0,0,0,1]
        imyy = rlimwt[:,0,0,0,1,1]
        wtxx = rlimwt[:,0,0,0,0,2]
        wtyy = rlimwt[:,0,0,0,1,2]
        rl = (rlxx*wtxx+rlyy*wtyy)/(wtxx+wtyy)
        im = (imxx*wtxx+imyy*wtyy)/(wtxx+wtyy)
        wt = 4.0/(1.0/wtxx + 1.0/wtyy)
        amp=np.sqrt(rl**2+im**2)

        deltaphase = 2*math.pi*(u*deltara + v*deltadec)
        phase = np.arctan(im/rl)+deltaphase
        rlcor = np.cos(phase)*amp
        rlcor[np.where( rl<0)] = -1.0*rlcor
        imcor = np.sin(phase)*amp
        imcor[np.where( rl<0)] = -1.0*imcor
    ampcor = np.sqrt(rlcor**2+imcor**2)


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
    for j in range(bins-1):
        if np.size(np.where( (depdist>=distrange[j]) & (depdist<=distrange[j+1]) ))==0:
            rlavg[j]=0
            imavg[j]=0
        else:
            rlavg[j] = np.sum((rlcor*wt)[np.where( (depdist>=distrange[j]) & (depdist<=distrange[j+1]) )]) / (np.sum(wt[np.where( (depdist>=distrange[j]) & (depdist<=distrange[j+1]) )]))

            imavg[j] = np.sum((imcor*wt)[np.where( (depdist>=distrange[j]) & (depdist<=distrange[j+1]) )]) / (np.sum(wt[np.where( (depdist>=distrange[j]) & (depdist<=distrange[j+1]) )]))
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
        modrlimwt = modfits[0].data['data']
        modheader = modfits[0].header
        modu = moddata['UU']
        modv = moddata['VV']
        modrlxx = np.ravel(modrlimwt[:,0,0,0,0,0], order='F')
        modrlyy = np.ravel(modrlimwt[:,0,0,0,1,0], order='F')
        modimxx = np.ravel(modrlimwt[:,0,0,0,0,1], order='F')
        modimyy = np.ravel(modrlimwt[:,0,0,0,1,1], order='F')
        modrl = (modrlxx*wtxx+modrlyy*wtyy)/(wtxx+wtyy)
        modim = (modimxx*wtxx+modimyy*wtyy)/(wtxx+wtyy)
        modfreq = modheader['CRVAL4']
        modamp = np.sqrt((modrl**2+modim**2))

        modu1 = (modu*np.cos(modpa)-modv*np.sin(modpa))*np.cos(modinc)
        modv1 = modu*np.sin(modpa)+modv*np.cos(modpa)
        modu2 = modu1*np.cos(-1*modpa)-modv1*np.sin(-1*modpa)
        modv2 = modu1*np.sin(-1*modpa)+modv1*np.cos(-1*modpa)
        moduu = modu2*modfreq*1e-03
        modvv = modv2*modfreq*1e-03
        moddepdist = (moduu**2 + modvv**2)**0.5

        if xlim:
            if binning == 'lin':
                moddistrange = np.linspace(np.min(moddepdist), xlim, bins)
            if binning == 'log':
                moddistrange = np.logspace(np.log10(np.min(moddepdist)), np.log10(xlim), bins)
        else:
            if binning == 'lin':
                moddistrange = np.linspace(np.min(moddepdist), np.max(moddepdist), bins)
            if binning == 'log':
                moddistrange = np.logspace(np.log10(np.min(moddepdist)), np.log10(np.max(moddepdist)), bins)
        modrlavg = np.zeros(bins)
        modimavg = np.zeros(bins)
        for i in range(bins-1):
            if np.size(np.wherdfits.writetoe( (moddepdfits.writetodist>=moddistrange[i]) & (moddepdist<=moddistrange[i+1]) ))==0:
                modrlavg[i]=0
                modimavg[i]=0
            else:
                modrlavg[i] = np.sum((modrl*wt)[np.where( (moddepdist>=moddistrange[i]) & (moddepdist<=moddistrange[i+1]) )]) / (np.sum(wt[np.where( (moddepdist>=moddistrange[i]) & (moddepdist<=moddistrange[i+1]) )]))

                modimavg[i] = np.sum((modim*wt)[np.where( (moddepdist>=moddistrange[i]) & (moddepdist<=moddistrange[i+1]) )]) / (np.sum(wt[np.where( (moddepdist>=moddistrange[i]) & (moddepdist<=moddistrange[i+1]) )]))


#Plotting
    plt.subplot(211)
    plt.errorbar( distrange, rlavg*1e03, yerr=rlSEM*1e03, marker='.', fmt=' ' )
    plt.axhline()
    if len(model)>4:
        plt.plot(moddistrange, modrlavg*1e03) 
    plt.title('Real Flux as a Function of Deprojected Distance')
    plt.xlabel('$k\lambda$')
    plt.ylabel('mJy')

    plt.subplot(212)
    plt.errorbar( distrange, imavg*1e03, yerr=imSEM*1e03, marker='.', fmt=' ' )
    plt.axhline()
    if len(model)>4:
        plt.plot(moddistrange, modimavg*1e03) 
    plt.title('Imaginary Flux as a Function of Deprojected Distance')
    plt.xlabel('$k\lambda$')
    plt.ylabel('mJy')
    #plt.savefig('49cetibinned.pdf')
    plt.show()

#Creating visibility fits file
if outputtype in ('im', 'both'):
    data['UU'] = u2
    data['VV'] = v2
    call('rm -rf '+image+'.*', shell=True)
    dfits.writeto(''+image+'.uvf')
    call('fits in='+image+'.uvf op=uvin out='+image+'.vis', shell=True)
    if datatype == 'cont':
        call('invert vis='+image+'.vis map='+image+'.mp beam='+image+'.bm cell=0.3 imsize=256 options=systemp,mfs robust=2', shell=True)
    if datatype == 'line':
        call('invert vis='+image+'.vis map='+image+'.mp beam='+image+'.bm cell=0.3 imsize=256 '+imchans+' robust=2', shell=True)
    call('clean map='+image+'.mp beam='+image+'.bm out='+image+'.cl niters=300', shell=True)
    call('restor map='+image+'.mp beam='+image+'.bm model='+image+'.cl out='+image+'.cm', shell=True)
    call('cgdisp in='+image+'.cm '+chandims+' device=/xs labtyp=arcsec beamtyp=b,l,4', shell=True)

