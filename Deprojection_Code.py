import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import math
import pyfits
from subprocess import call
import os
plt.close()

#INPUT VARIABLES
print "This code will read in a FITS file and deproject its visibilities. It can handle line emission and continuum data, and can radially average continuum to create a plot of flux as a function of radially deprojected distance. It can also image the disk to create either a single continuum image, a zeroth moment map, or channel maps."
print
fileinp = raw_input('Please enter the filename: ')
print
print "If there is an offset between the center of the source and the center of the image, please enter the shift required to center the source. (if the source center is at (-1.0,-1.0) you would enter (1.0,1.0)"
print
deltara = math.radians( input('RA shift in arcsecs: ')/3600.0)
deltadec = math.radians( input('Dec shift in arcsecs: ')/3600.0) 
print
pa = math.radians( input('Please enter the SKY PA in degrees: '))
inc = math.radians( input('Please enter the inclination of the disk in degrees: '))
print
outputtype = raw_input('Please enter \'rad\' for radial averaging or \'im\' for imaging: ')
if outputtype == 'im':
   image = raw_input('Please enter the filename you would like your image(s) to have, not including the filetype: ')
if outputtype == 'rad':
    bins = input('Please enter the number of bins you would like: ')
    binning = raw_input('Please enter \'log\' for logarithmic binning or \'lin\' for linear binning: ')
    print
    try:
       xlim = input('If you would like to set a value for the maximum distance plotted, please enter it in $k\lambda$ here. Otherwise, press return: ')
    except SyntaxError:
       xlim = None
    try:
        print
        savefig = raw_input('If you would like to save the radially averaged plot, please enter the filename (including extension) that you would like it to have. Otherwise, press return: ')
    except SyntaxError:
        savefig = None
dfits = fits.open(fileinp)
header = dfits[0].header
if header['NAXIS4'] == 1:
    datatype = 'cont'
else:
    datatype = 'line'
if datatype == 'line':
    if outputtype == 'rad':
        momchans = 'moment'
    else:
        print
        momchans = raw_input('Please enter \'moment\' if you would like to create a moment map, or \'channel\' if you would like to create channel maps: ')
        print
    if momchans == 'channel':
        print
        chandims = 'nxy=' + raw_input('Please enter the dimensions of the channel maps to be shown on each page, in the form \'columns,rows\': ')
if outputtype == 'rad':
    try:   
        print
        model = raw_input('If you would like to compare with a model, please include its filename here. Otherwise, press return: ')
    except SyntaxError:
       model = None
#fileinp='ceti250ch.uvf'
#outputtype='rad'
#bins=15
#binning='log'
#spa=108.5
#pa=math.radians(spa)
#inc=math.radians(79.1)
#model='vf'
#image='49cetichannelmap'
#deltara=0
#deltadec=0
#xlim = 300
#datatype= 'line'
#chansin = 'line=channel,'+'23,1,1,1'
#chandims = 'nxy=6,3'
#chans = range(0,2)
#momchans = 'aaaaaa'    


#Reading in fits file/image (assuming ALMA data)
data = dfits[0].data
rlimwt = dfits[0].data['data']
chans = range(header['NAXIS4'])
if datatype == 'cont':
    rlxx = rlimwt[:,0,0,:,0,0]
    rlyy = rlimwt[:,0,0,:,1,0]
    imxx = rlimwt[:,0,0,:,0,1]
    imyy = rlimwt[:,0,0,:,1,1]
    wtxx = rlimwt[:,0,0,:,0,2]
    wtyy = rlimwt[:,0,0,:,1,2]
if datatype == 'line':
    rlxx = rlimwt[:,0,0,0,:,0,0]
    rlyy = rlimwt[:,0,0,0,:,1,0]
    imxx = rlimwt[:,0,0,0,:,0,1]
    imyy = rlimwt[:,0,0,0,:,1,1]
    wtxx = rlimwt[:,0,0,0,:,0,2]
    wtyy = rlimwt[:,0,0,0,:,1,2]
amp_xx=np.sqrt(rlxx**2+imxx**2)
amp_yy=np.sqrt(rlyy**2+imyy**2)
freq = header['CRVAL4']
u = data['UU']
v = data['VV']
visnum=len(u)


#Calculating phase shifts
deltaphase = 2*math.pi*(u*deltara + v*deltadec)
deltaphase = (np.tile(deltaphase, (len(chans), 1))).T
phase_xx = np.arctan(imxx/rlxx)+deltaphase
phase_yy = np.arctan(imyy/rlyy)+deltaphase
rlcor_xx = np.cos(phase_xx)*amp_xx
rlcor_yy = np.cos(phase_yy)*amp_yy
rlcor_xx[np.where( rlxx<0 )] = -1.0*rlcor_xx[np.where( rlxx<0 )]
rlcor_yy[np.where( rlyy<0 )] = -1.0*rlcor_yy[np.where( rlyy<0 )]
imcor_xx = np.cos(phase_xx)*amp_xx
imcor_yy = np.cos(phase_yy)*amp_yy
imcor_xx[np.where( imxx<0 )] = -1.0*imcor_xx[np.where( imxx<0 )]
imcor_yy[np.where( imyy<0 )] = -1.0*imcor_yy[np.where( imyy<0 )]

rlcor = np.array((rlcor_xx*wtxx+rlcor_yy*wtyy)/(wtxx+wtyy))
imcor = np.array((imcor_xx*wtxx+imcor_yy*wtyy)/(wtxx+wtyy))
wtcor = 4.0/(1.0/wtxx + 1.0/wtyy)


#Rotation of (u,v) coordinates back to PA=0 and deprojection
u1 = (u*np.cos(pa)-v*np.sin(pa))*np.cos(inc)
v1 = (u*np.sin(pa)+v*np.cos(pa))
#Rotation of (u,v) coordinates back to sky PA 
u2 = u1*np.cos(-1.0*pa)-v1*np.sin(-1.0*pa)
v2 = u1*np.sin(-1.0*pa)+v1*np.cos(-1.0*pa)
uu = u2*freq*1e-03
vv = v2*freq*1e-03
depdist = (uu**2 + vv**2)**0.5


#Radial Averaging
if outputtype == 'rad':
    if datatype == 'line':
        for i in range(visnum):
            rlcor[i] = np.sum(rlcor[i,:])
            imcor[i] = np.sum(imcor[i,:])
            wtcor[i] = np.sum(wtcor[i,:])
    if xlim:
        if binning == 'lin':
            distrange = np.linspace(np.min(depdist), xlim, bins)
        elif binning == 'log':
            distrange = np.logspace(np.log10(np.min(depdist)), np.log10(xlim), bins)
    else:
        if binning == 'lin':   
            distrange = np.linspace(np.min(depdist), np.max(depdist), bins)
        elif binning == 'log':
            distrange = np.logspace(np.log10(np.min(depdist)), np.log10(np.max(depdist)), bins)
    rlavg = np.zeros(bins)
    imavg = np.zeros(bins)
    for i in range(bins-1):
        if np.size(np.where( (depdist>=distrange[i]) & (depdist<=distrange[i+1]) ))==0:
            rlavg[i]=0
            imavg[i]=0
        else:
            rlavg[i] = np.sum((rlcor*wtcor)[np.where( (depdist>=distrange[i]) & (depdist<=distrange[i+1]) )]) / (np.sum(wtcor[np.where( (depdist>=distrange[i]) & (depdist<=distrange[i+1]) )]))

            imavg[i] = np.sum((imcor*wtcor)[np.where( (depdist>=distrange[i]) & (depdist<=distrange[i+1]) )]) / (np.sum(wtcor[np.where( (depdist>=distrange[i]) & (depdist<=distrange[i+1]) )]))
    ampavg = np.sqrt(rlavg**2 + imavg**2)


    #Standard Error of the Mean
    rlsamplestandev = np.zeros(bins)
    imsamplestandev = np.zeros(bins)
    N = np.zeros(bins)
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
        modchans = range(modheader['NAXIS4'])
        if datatype == 'cont':
            modrl = modrlimwt[::2,0,0,:,0,0]
            modim = modrlimwt[::2,0,0,:,0,1]
        elif datatype == 'line':
            modrl = modrlimwt[::2,0,0,0,:,0,0]
            modim = modrlimwt[::2,0,0,0,:,0,1]
        modfreq = modheader['CRVAL4']
        modu = moddata['UU'][::2]
        modv = moddata['VV'][::2]
        modvisnum=len(modu)
        modfreq = modheader['CRVAL4']

        modu1 = (modu*np.cos(pa)-modv*np.sin(pa))*np.cos(inc)
        modv1 = modu*np.sin(pa)+modv*np.cos(pa)
        modu2 = modu1*np.cos(-1.0*pa)-modv1*np.sin(-1.0*pa)
        modv2 = modu1*np.sin(-1.0*pa)+modv1*np.cos(-1.0*pa)
        moduu = modu2*modfreq*1e-03
        modvv = modv2*modfreq*1e-03
        moddepdist = (moduu**2 + modvv**2)**0.5

        if xlim:
            if binning == 'lin':
                moddistrange = np.linspace(np.min(moddepdist), xlim, bins)
            elif binning == 'log':
                moddistrange = np.logspace(np.log10(np.min(moddepdist)), np.log10(xlim), bins)
        else:
            if binning == 'lin':
                moddistrange = np.linspace(np.min(moddepdist), np.max(moddepdist), bins)
            elif binning == 'log':
                moddistrange = np.logspace(np.log10(np.min(moddepdist)), np.log10(np.max(moddepdist)), bins)
        modrlavg = np.zeros(bins)
        modimavg = np.zeros(bins)
        for i in range(bins-1):
            if np.size(np.where( (moddepdist>=moddistrange[i]) & (moddepdist<=moddistrange[i+1]) ))==0:
                modrlavg[i]=0
                modimavg[i]=0
            else:
                modrlavg[i] = np.sum((modrl*wtcor)[np.where( (moddepdist>=moddistrange[i]) & (moddepdist<=moddistrange[i+1]) )]) / (np.sum(wtcor[np.where( (moddepdist>=moddistrange[i]) & (moddepdist<=moddistrange[i+1]) )]))

                modimavg[i] = np.sum((modim*wtcor)[np.where( (moddepdist>=moddistrange[i]) & (moddepdist<=moddistrange[i+1]) )]) / (np.sum(wtcor[np.where( (moddepdist>=moddistrange[i]) & (moddepdist<=moddistrange[i+1]) )]))


#Plotting
    plt.subplot(211)
    plt.errorbar( distrange, rlavg*1e03, yerr=rlSEM*1e03, marker='.', fmt=' ' )
    plt.axhline()
    if len(model)>4:
        plt.plot(moddistrange, modrlavg*1e03) 
    plt.title('Flux as a Function of Deprojected Distance')
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
    if savefig:
        plt.savefig(''+savefig+'')
    plt.show()


#Creating visibility fits file
if outputtype == 'im':
    data['UU'] = u2
    data['VV'] = v2

    call('rm -rf '+image+'.*', shell=True)
    dfits.writeto(''+image+'.uvf')
    call('fits in='+image+'.uvf op=uvin out='+image+'.vis', shell=True)
    if datatype == 'cont':
        call('invert vis='+image+'.vis map='+image+'.mp beam='+image+'.bm cell=0.3 imsize=256 options=systemp,mfs robust=2', shell=True)
    elif datatype == 'line':
            call('invert vis='+image+'.vis map='+image+'.mp beam='+image+'.bm cell=0.3 imsize=256 robust=2', shell=True)
    call('clean map='+image+'.mp beam='+image+'.bm out='+image+'.cl niters=500', shell=True)
    call('restor map='+image+'.mp beam='+image+'.bm model='+image+'.cl out='+image+'.cm', shell=True)
    if datatype == 'cont':
        call('cgdisp in='+image+'.cm device=/xs labtyp=arcsec beamtyp=b,l,4', shell=True)
    elif datatype == 'line':
        if momchans == 'moment':
            call('moment in='+image+'.cm mom=0 out='+image+'.m0', shell=True)
            call('cgdisp in='+image+'.m0 device=/xs labtyp=arcsec beamtyp=b,l,4', shell=True)
        else:  
            call('cgdisp in='+image+'.cm '+chandims+' device=/xs labtyp=arcsec beamtyp=b,l,4', shell=True)

    

