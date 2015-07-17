from matplotlib.ticker import AutoMinorLocator,LinearLocator,NullLocator
import matplotlib.gridspec as gridspec
import math
import numpy as np
import matplotlib.pyplot as plt
import pyfits
import pylab
from matplotlib.patches import Ellipse
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patheffects as PathEffects

#.fits inputs
dat = ['49Ceti_.75arcsectaper.fits','49Ceti_robust2.fits']

#make the colormap
def make_cmap(colors, position=None, bit=False):
	import matplotlib as mpl
	bit_rgb = np.linspace(0,1,256)
	if position == None:
		position = np.linspace(0,1,len(colors))
	else:
		if len(position) != len(colors):
			sys.exit("position length must be the same as colors")
		elif position[0] != 0 or position[-1] != 1:
			sys.exit('position must start with 0 and end with 1')
	if bit:
		for i in range(len(colors)):
			colors[i] = (bit_rgb[colors[i][0]],					                             bit_rgb[colors[i][1]],					                             bit_rgb[colors[i][2]])
	cdict = {'red':[], 'green':[], 'blue':[]}
	for pos, color in zip(position, colors):
		cdict['red'].append((pos, color[0], color[0]))
		cdict['green'].append((pos, color[1], color[1]))
		cdict['blue'].append((pos, color[2], color[2]))
	cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
	return cmap

colors = [(255,255,255),(255,255,255),(248,255,255),(240,255,255),(210,253,255),(184,252,255),(192,244,204),(155,255,145),(210,200,12),(230,180,7),(236,124,13),(233,100,25),(230,60,30),(228,30,45),(227,7,63),(222,5,150),(218,2,218)]# FINAL!

my_cmap = make_cmap(colors,bit=True)

#set image size: 7.5x4.2 makes square frames with the colorbar at its current size
fig = plt.figure(figsize=(7.5,4.2))
outer_grid = gridspec.GridSpec(1,2,wspace = 0.0, hspace = 0.0)
pylab.xticks([])
pylab.yticks([])
import matplotlib.cm as cm

def plotcmd(dat,ax,levs1):
	hdulist = pyfits.open(dat[i])
	imagedat = hdulist[0].data
	imagedatsqueeze = imagedat.squeeze()
	imgGood = imagedatsqueeze[:-1,:]
	cell=.009765625
	col = ax.pcolormesh(imgGood,cmap=my_cmap,vmin=np.min(imgGood),vmax=np.max(imgGood))
	#col = ax.pcolormesh(imgGood,cmap=my_cmap,vmin=0,vmax=np.max(imgGood)) # you can scale from the zero of the map instead of the min if you want-- with the current colormap, scaling from the minimum works better
	ax.contour(imgGood, levels=levs1, linewidths=.75,colors='k')
	ax.contour(imgGood, -levs1, linewidths=.75,colors='k')
	ax.xaxis.set_major_locator(LinearLocator(11))
	ax.yaxis.set_major_locator(LinearLocator(11))
	ax.minorticks_on()
	ax.tick_params(which='minor',length=2, width=1, colors='k')
	ax.tick_params(which='major',length=5, width=1, colors='k')
	ax.tick_params(axis='x',labelcolor='k')
	ax.tick_params(axis='y',labelcolor='k')
	ax.xaxis.label.set_color('k')
	ax.yaxis.label.set_color('k')
	ax.set_xticklabels([])
	ax.set_yticklabels([])
	#describe the beams for each image
	if i == 0:
		beam = Ellipse(xy=(80,70),height=.54/.009765625,width=0.43/.009765625,angle=-8.49e1,facecolor='none',edgecolor='k',hatch='//')
	if i == 1:
		beam = Ellipse(xy=(90,80),height=.91/.009765625,width=0.77/.009765625,angle=80.1,facecolor='none',edgecolor='k',hatch='//')
	ax.add_patch(beam)
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("top", size="8%", pad=0.0)
	colbar = fig.colorbar(col,cax=cax,orientation='horizontal')
	axT = colbar.ax.twiny()
	print np.min(imgGood)
	print np.max(imgGood)
	axT.set_xlim(np.min(imgGood),np.max(imgGood))
	#set x colbar ticks and labels
	if i == 0:
		axT.set_xticks(np.arange(-.0002,0.0007001,0.0002))
		axT.set_xticklabels(range(-200,701,200),fontsize=7,rotation=30)
		axT.tick_params(axis='both',which='major',pad=1)
	if i == 1:
		axT.set_xticks(np.arange(-.0002,0.0014001,0.0002))
		axT.set_xticklabels(range(-200,1401,200),fontsize=7,rotation=30)
		axT.tick_params(axis='both',which='major',pad=1)
	axT.set_yticks([])
	axT.set_yticklabels('',visible=False)
	colbar.set_ticks([])
	colbar.set_ticklabels([])
	colbar.ax.set_xticklabels('',visible=False)
	if i == 1:
		ax.text(705,967,'0.75" gaussian',fontsize=8,path_effects=[PathEffects.withStroke(linewidth=2,foreground="w")])
		ax.text(825,920,'taper',fontsize=8,path_effects=[PathEffects.withStroke(linewidth=2,foreground="w")])
	if ax.is_first_col():
		ax.set_xlabel(r'$\Delta \alpha$ (")',fontsize=10,labelpad=1)
		ax.set_ylabel(r'$\Delta \delta$ (")',fontsize=10,labelpad=1)
		x = ['','-4','','-2','','0','','2','','4','']
		ax.set_xticklabels(x)
		ax.set_yticklabels(x)
		for tick in ax.xaxis.get_major_ticks():
			tick.label.set_fontsize(8)
		for tick in ax.yaxis.get_major_ticks():
			tick.label.set_fontsize(8)
		line = pylab.Line2D((615,785.666),(100,100),lw=1,color='k') ##100AU scale bar
		ax.add_line(line)
		line = pylab.Line2D((615,615),(95,105),lw=1,color='k')
		ax.add_line(line)
		line = pylab.Line2D((785.666,785.666),(95,105),lw=1,color='k')
		ax.add_line(line)
		ax.text(800,967,'robust=2',fontsize=8,path_effects=[PathEffects.withStroke(linewidth=2,foreground="w")])
		ax.text(627,115,'100 AU',fontsize=8,path_effects=[PathEffects.withStroke(linewidth=2,foreground="w")])
		ax.text(444,967,'49 Ceti',fontsize=8,path_effects=[PathEffects.withStroke(linewidth=2,foreground="w")])#color='k')
		ax.text(381,920,r'ALMA 850$\mu m$',fontsize=8,path_effects=[PathEffects.withStroke(linewidth=2,foreground="w")])
		ax.text(444,1053,r'$\mu Jy / bm$',fontsize=8,path_effects=[PathEffects.withStroke(linewidth=2,foreground="w")])


for i in range(2):
	if i == 0:
		levels = np.array([j*5.55e-5 for j in range(2,40,2)])
	if i == 1:
		levels = np.array([j*7.50e-5 for j in range(2,50,2)])
	inner_grid = gridspec.GridSpecFromSubplotSpec(1,2, subplot_spec = outer_grid[i], wspace = 0.0, hspace = 0.0)
	ax = plt.Subplot(fig, outer_grid[i])
	plotcmd(dat,ax, levels)
	fig.add_subplot(ax)

pylab.gcf().subplots_adjust(bottom=0.15)
#plt.savefig("/home/jliemansifry/Desktop/49CET_doubleplot_FINALCOLORS.png",dpi=400)
plt.show()
