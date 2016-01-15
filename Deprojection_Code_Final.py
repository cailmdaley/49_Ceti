from astropy.io import fits
import os
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
plt.close()

#INPUT VARIABLES//Opening fits file
print
print "This code will read in a FITS file (either data or a model//residuals) and deproject its visibilities. It can handle line emission and continuum data, and can average data to create a plot of flux as a function of radially deprojected distance. It can also image the disk to create either a single continuum image, a zeroth moment map, or channel maps."
print
try:
    fileinp = raw_input('If you would like to deproject observed data, please enter its filename here. Otherwise, please hit return: ')
except SyntaxError:
    fileinp = None
try:
    print
    model = raw_input('If you would like to deproject a model//residuals, please include its filename here. Otherwise, press return. (Note: If you attempt to image both observed data and a model, only the model will be imaged) ')
except SyntaxError:
    model = None
    #Reading in fits file//determining channels
if fileinp:
    dfits = fits.open(fileinp)
    if dfits[0].header['NAXIS4'] == 1:
        datatype = 'cont'
    else:
        datatype = 'line'
if model:
        modfits = fits.open(model)
        if modfits[0].header['NAXIS4'] == 1:
            datatype = 'cont'
        else:
            datatype = 'line'
if fileinp:
    print
    print "If there is an offset between the center of the source and the center of the image, please enter the shift required to center the source. (if the source center is at (-1.0,-1.0) you would enter (1.0,1.0)"
    print
    deltara = math.radians( input('RA shift in arcsecs: ')/3600.0)
    deltadec = math.radians( input('Dec shift in arcsecs: ')/3600.0)
print
pa = math.radians(-71.1505145) #input('Please enter the SKY PA in degrees: '))
inc = math.radians( 79.2284311) #input('Please enter the inclination of the disk in degrees: '))
print
outputtype = raw_input('Please enter \'rad\' for radial averaging or \'im\' for imaging: ')
if outputtype == 'im':
    image = raw_input('Please enter the filename you would like your image(s) to have, not including the extension: ')
    imsize = input('Please how many arcseconds across you would like your image to be: ')
    pixsize = str(imsize/1024.0)
    if datatype =='cont':
        cutoff = raw_input('Please enter the flux level at which you would like to stop cleaning, in Jy: ')
if outputtype == 'rad':
    binning = raw_input('Please enter \'log\' for logarithmic binning or \'lin\' for linear binning: ')
    bins = input('Please enter the number of bins you would like: ')
    print
    try:
       xmax = input('If you would like to set a value for the maximum distance plotted, please enter it in $k\lambda$ here. Otherwise, press return: ')
    except SyntaxError:
       xmax = None
    try:
        print
        savefig = raw_input('If you would like to save the radially averaged plot, please enter the filename (excluding extension) that you would like it to have. Otherwise, press return: ')
    except SyntaxError:
        savefig = None
if datatype == 'line':
    if outputtype == 'im':
        print
        momchans = raw_input('Please enter \'moment\' if you would like to create a moment map, or \'channel\' if you would like to create channel maps: ')
        if momchans == 'channel':
            imchans = 'line=channel,' + raw_input ( 'Please enter the channels you would like to image in standard line MIRIAD format (\'number of channels,starting channel index,channel width,step size\'). Indexing starts at 1: ')
            chandims = 'nxy=' + raw_input('Please enter the dimensions of the displayed channel maps in the form of \'x,y\': ')

#Reading in polarized real & imaginary visibilities & weights (assuming ALMA data)
if fileinp:
    rlimwt = dfits[0].data['data']
    chans = range(dfits[0].header['NAXIS4'])
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
    u = dfits[0].data['UU']
    v = dfits[0].data['VV']
    visnum=len(u)

    #Calculating & correcting for phase shifts
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

    #Taking weighted average of polarized visibilities & weights
    rlcor = np.array(0.25*(1.0/wtxx + 1.0/wtyy)*(rlcor_xx*wtxx + rlcor_yy*wtyy))
    imcor = np.array(0.25*(1.0/wtxx + 1.0/wtyy)*(imcor_xx*wtxx + imcor_yy*wtyy))
    wtcor = 4.0/(1.0/wtxx + 1.0/wtyy)

    #Rotation of (u,v) coordinates back to PA=0 and deprojection
    u1 = (u*np.cos(pa)-v*np.sin(pa))*np.cos(inc)
    v1 = (u*np.sin(pa)+v*np.cos(pa))

    #Rotation of (u,v) coordinates back to sky PA
    u2 = u1*np.cos(-1.0*pa)-v1*np.sin(-1.0*pa)
    v2 = u1*np.sin(-1.0*pa)+v1*np.cos(-1.0*pa)

    #Converting to kilolambda and Jy*km/s
    freq = dfits[0].header['CRVAL4']
    dfreq= dfits[0].header['CDELT4']
    chanwidth = 2.99792458e05*-1.0*dfreq/freq
    uu = u2*freq*1e-03
    vv = v2*freq*1e-03
    if outputtype == 'rad':
        if datatype == 'line':
            for i in range(visnum):
                rlcor[i] = np.sum(rlcor[i,:])*chanwidth
                imcor[i] = np.sum(imcor[i,:])*chanwidth
                wtcor[i] = np.sum(wtcor[i,:])*chanwidth
            rlcor = np.delete(rlcor, np.s_[1:], 1)
            imcor = np.delete(imcor, np.s_[1:], 1)
            wtcor = np.delete(wtcor, np.s_[1:], 1)

    #Radial Averaging
        depdist = (uu**2 + vv**2)**0.5
        if xmax:
            if binning == 'lin':
                binbounds = np.linspace(np.min(depdist), xmax, bins+1)
            elif binning == 'log':
                binbounds = np.logspace(np.log10(np.min(depdist)), np.log10(xmax), bins+1)
        else:
            xmax = max(depdist)
            if binning == 'lin':
                binbounds = np.linspace(np.min(depdist), np.max(depdist), bins+1)
            elif binning == 'log':
                binbounds = np.logspace(np.log10(np.min(depdist)), np.log10(np.max(depdist)), bins+1)
        binctr = np.zeros(bins)
        rlavg = np.zeros(bins)
        imavg = np.zeros(bins)
        rlwtedSEM = np.zeros(bins)
        imwtedSEM = np.zeros(bins)
        N = np.zeros(bins)
        rlxbar = np.zeros(bins)
        rlterm1 = np.zeros(bins)
        rlterm2 = np.zeros(bins)
        rlterm3 = np.zeros(bins)
        imxbar = np.zeros(bins)
        imterm1 = np.zeros(bins)
        imterm2 = np.zeros(bins)
        imterm3 = np.zeros(bins)
        for i in range(bins):
            binctr[i] = (binbounds[i] + binbounds[i+1])/2.0
            if datatype == 'cont':
                rlavg[i] = (np.sum((rlcor*wtcor)[np.where( (depdist>=binbounds[i]) & (depdist<=binbounds[i+1]) )]) / (np.sum(wtcor[np.where( (depdist>=binbounds[i]) & (depdist<=binbounds[i+1]) )])))*1e03
                imavg[i] = (np.sum((imcor*wtcor)[np.where( (depdist>=binbounds[i]) & (depdist<=binbounds[i+1]) )]) / (np.sum(wtcor[np.where( (depdist>=binbounds[i]) & (depdist<=binbounds[i+1]) )])))*1e03
            else:
                rlavg[i] = np.sum((rlcor*wtcor)[np.where( (depdist>=binbounds[i]) & (depdist<=binbounds[i+1]) )]) / (np.sum(wtcor[np.where( (depdist>=binbounds[i]) & (depdist<=binbounds[i+1]) )]))
                imavg[i] = np.sum((imcor*wtcor)[np.where( (depdist>=binbounds[i]) & (depdist<=binbounds[i+1]) )]) / (np.sum(wtcor[np.where( (depdist>=binbounds[i]) & (depdist<=binbounds[i+1]) )]))


            #Standard Error of the Mean
            N[i] = np.size(np.where( (depdist>=binbounds[i]) & (depdist<=binbounds[i+1]) ))
            rlxbar[i] = np.sum((rlcor*wtcor)[np.where((depdist>=binbounds[i]) & (depdist<=binbounds[i+1]))])/np.sum(wtcor[np.where((depdist>=binbounds[i]) & (depdist<=binbounds[i+1]))])
            imxbar[i] = np.sum((imcor*wtcor)[np.where((depdist>=binbounds[i]) & (depdist<=binbounds[i+1]))])/np.sum(wtcor[np.where((depdist>=binbounds[i]) & (depdist<=binbounds[i+1]))])
            rlterm1[i] = np.sum(((rlcor*wtcor)[np.where((depdist>=binbounds[i]) & (depdist<=binbounds[i+1]))]-np.mean(wtcor[np.where((depdist>=binbounds[i]) & (depdist<=binbounds[i+1]))])*rlxbar[i])**2)
            rlterm2[i] = -2*rlxbar[i]*np.sum((wtcor[np.where((depdist>=binbounds[i]) & (depdist<=binbounds[i+1]))]-np.mean(wtcor[np.where((depdist>=binbounds[i]) & (depdist<=binbounds[i+1]))]))*((wtcor*rlcor)[np.where((depdist>=binbounds[i]) & (depdist<=binbounds[i+1]))]-np.mean(wtcor[np.where((depdist>=binbounds[i]) & (depdist<=binbounds[i+1]))])*rlxbar[i]))
            rlterm3[i] = rlxbar[i]**2*np.sum((wtcor[np.where((depdist>=binbounds[i]) & (depdist<=binbounds[i+1]))]-np.mean(wtcor[np.where((depdist>=binbounds[i]) & (depdist<=binbounds[i+1]))]))**2)
            imterm1[i] = np.sum(((imcor*wtcor)[np.where((depdist>=binbounds[i]) & (depdist<=binbounds[i+1]))]-np.mean(wtcor[np.where((depdist>=binbounds[i]) & (depdist<=binbounds[i+1]))])*imxbar[i])**2)
            imterm2[i] = -2*imxbar[i]*np.sum((wtcor[np.where((depdist>=binbounds[i]) & (depdist<=binbounds[i+1]))]-np.mean(wtcor[np.where((depdist>=binbounds[i]) & (depdist<=binbounds[i+1]))]))*((wtcor*imcor)[np.where((depdist>=binbounds[i]) & (depdist<=binbounds[i+1]))]-np.mean(wtcor[np.where((depdist>=binbounds[i]) & (depdist<=binbounds[i+1]))])*imxbar[i]))
            imterm3[i] = imxbar[i]**2*np.sum((wtcor[np.where((depdist>=binbounds[i]) & (depdist<=binbounds[i+1]))]-np.mean(wtcor[np.where((depdist>=binbounds[i]) & (depdist<=binbounds[i+1]))]))**2)
            rlwtedSEM[i] = np.sqrt((N[i]/((N[i]-1)*np.sum(wtcor[np.where((depdist>=binbounds[i]) & (depdist<=binbounds[i+1]))])**2))*(rlterm1[i]+rlterm2[i]+rlterm3[i]))
            imwtedSEM[i] = np.sqrt((N[i]/((N[i]-1)*np.sum(wtcor[np.where((depdist>=binbounds[i]) & (depdist<=binbounds[i+1]))])**2))*(imterm1[i]+imterm2[i]+imterm3[i]))
        if datatype == 'cont':
            rlwtedSEM = rlwtedSEM*1e03
            imstedSEM = imwtedSEM*1e03


#Reading in model visibilities in Stokes I
if model:
    if datatype == 'cont':
        modrl = modfits[0].data['data'][:,0,0,:,0,0].ravel()*1e03
        modim = modfits[0].data['data'][:,0,0,:,0,1].ravel()*1e03
    elif datatype == 'line':
        modrl = modfits[0].data['data'][:,0,0,0,:,0,0].ravel()
        modim = modfits[0].data['data'][:,0,0,0,:,0,1].ravel()
    modu = modfits[0].data['UU']
    modv = modfits[0].data['VV']
    modvisnum=len(modu)
    modfreq = modfits[0].header['CRVAL4']

    #Deprojecting model visibilities
    modu1 = (modu*np.cos(pa)-modv*np.sin(pa))*np.cos(inc)
    modv1 = modu*np.sin(pa)+modv*np.cos(pa)
    modu2 = modu1*np.cos(-1.0*pa)-modv1*np.sin(-1.0*pa)
    modv2 = modu1*np.sin(-1.0*pa)+modv1*np.cos(-1.0*pa)

    #Converting to kilolambda and Jy*km/s
    moduu = modu2*modfreq*1e-03
    modvv = modv2*modfreq*1e-03
    if outputtype == 'rad':
        if datatype == 'line':
            for i in range(visnum):
                modrl[i] = np.sum(modrl[i,:])*chanwidth
                modim[i] = np.sum(modim[i,:])*chanwidth
            modrl= np.delete(modrl, np.s_[1:], 1)
            modim= np.delete(modim, np.s_[1:], 1)

    #Sorting model visibilities by uv distance
    moddepdist = (moduu**2 + modvv**2)**0.5
    moddistinds = moddepdist.argsort()
    moddepdist = moddepdist[moddistinds]
    modrl = modrl[moddistinds]


#Plotting: setting y scales to be the same for both plots
if outputtype == 'rad':
    if fileinp:
        rlmax = np.ceil(max(rlavg) + rlwtedSEM[np.where(rlavg==max(rlavg))])
        rlmin = np.floor(min(rlavg) - rlwtedSEM[np.where(rlavg==min(rlavg))])
        imrange = np.ceil(max(np.absolute(imavg)) + max(imwtedSEM))
    if model:
        if fileinp:
            if max(modrl) >= rlmax:
                rlmax = np.ceil(max(modrl))
            if min(modrl) <= rlmin:
                rlmin = np.floor(min(modrl))
        else:
            rlmax = np.ceil(max(modrl))
            rlmin = np.floor(min(modrl))
    rlrange = rlmax - rlmin

#Plotting: setting up gridspec, setting ticks and labels, etc.
    if fileinp:
        gs = gridspec.GridSpec(2, 1, height_ratios=[rlrange,2*imrange])
    else:
        gs = gridspec.GridSpec(1, 1)

    plt.subplot(gs[0])
    if model:
        plt.plot(moddepdist, modrl, linewidth=2, c='g')
        plt.xticks(np.arange(0, max(moddepdist)+.0001, 20))
    plt.xlim(0, xmax)
    plt.axhline(c='k', ls='--', linewidth=2)
    if datatype == 'cont':
        plt.ylabel('Re (mJy)', fontsize=15)
    else:
        plt.ylabel('Re (Jy$\mathcal{*}$km/s)', fontsize=15)
    if fileinp:
        plt.errorbar( binctr, rlavg, yerr=rlwtedSEM, marker='.', c='b', fmt=' ' )
        plt.xticks(np.arange(0, max(binbounds)+.0001, 20))

        plt.subplot(gs[1])
        plt.errorbar( binctr, imavg, yerr=imwtedSEM, marker='.', fmt=' ' )
        plt.xlim(0, xmax)
        plt.axhline(c='k', ls='--', linewidth=2)
        plt.xticks(np.arange(0, max(binbounds)+.0001, 20))
        if datatype == 'cont':
            plt.yticks(np.linspace(-1*imrange, imrange, 3))
        else:
            plt.yticks(np.linspace(np.floor((min(imavg) - imwtedSEM[np.where(imavg==min(imavg))])), np.ceil((max(imavg) + imwtedSEM[np.where(imavg==max(imavg))])), 3))
        plt.xlabel(r'$\mathcal{R}_{uv}$ ($k\lambda$)', fontsize=15, labelpad=10)
        if datatype == 'cont':
            plt.ylabel('Im (mJy)', fontsize=15)
        else:
            plt.ylabel(r'Im (Jy$\mathcal{*}$km/s)', fontsize=15)
    if savefig:
        plt.savefig(''+savefig+'.eps')
    plt.show()

#Reading data back into fits file, cleaning in MIRIAD, converting back to fits file
if outputtype == 'im':
    if fileinp:
        dfits[0].data['UU'] = u2
        dfits[0].data['VV'] = v2
        if datatype == 'cont':
            rlimwt[:,0,0,0,0,0]=np.ravel(rlcor_xx)
            rlimwt[:,0,0,0,1,0]=np.ravel(rlcor_yy)
            rlimwt[:,0,0,0,0,1]=np.ravel(imcor_xx)
            rlimwt[:,0,0,0,1,1]=np.ravel(imcor_yy)
        if datatype == 'line':
            rlimwt[:,0,0,0,:,0,0]=rlcor_xx
            rlimwt[:,0,0,0,:,1,0]=rlcor_yy
            rlimwt[:,0,0,0,:,0,1]=imcor_xx
            rlimwt[:,0,0,0,:,1,1]=imcor_yy
        dfits[0].data['data']=rlimwt
        os.system('rm -rf '+image+'.*')
        dfits.writeto(''+image+'.uvf')

    if model:
        modfits[0].data['UU'] = modu2
        modfits[0].data['VV'] = modv2
        os.system('rm -rf '+image+'.*')
        modfits.writeto(''+image+'.uvf')

    os.system('fits in='+image+'.uvf op=uvin out='+image+'.vis')
    if datatype == 'cont':
        os.system('invert vis='+image+'.vis map='+image+'.mp beam='+image+'.bm cell='+pixsize+' imsize=1024 options=systemp,mfs robust=2')
    elif datatype == 'line':
        if momchans == 'moment':
            os.system('invert vis='+image+'.vis map='+image+'.mp beam='+image+'.bm cell='+pixsize+' imsize=1024 robust=2')
        elif momchans == 'channel':
            os.system('invert vis='+image+'.vis map='+image+'.mp beam='+image+'.bm cell='+pixsize+' '+imchans+' imsize=1024 robust=2')
    if datatype == 'cont':
        os.system('clean map='+image+'.mp beam='+image+'.bm out='+image+'.cl cutoff='+cutoff+' niters=10000')
    elif datatype == 'line':
        os.system('clean map='+image+'.mp beam='+image+'.bm out='+image+'.cl niters=100')
    os.system('restor map='+image+'.mp beam='+image+'.bm model='+image+'.cl out='+image+'.cm')
    if datatype == 'cont':
            os.system('cgdisp in='+image+'.cm device=/xs labtyp=arcsec beamtyp=b,l,4')
    elif datatype == 'line':
        if momchans == 'moment':
            os.system('moment in='+image+'.cm mom=0 out='+image+'.m0')
            os.system('cgdisp in='+image+'.m0 device=/xs labtyp=arcsec beamtyp=b,l,4')
        else:
            os.system('cgdisp in='+image+'.cm '+chandims+' device=/xs labtyp=arcsec beamtyp=b,l,4')

