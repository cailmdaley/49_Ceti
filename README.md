# 49-Ceti
49 Ceti, a 40 million year old A1V type star 61 parsecs from Earth, is one of only a handful
of systems known to retain its molecular gas well into the debris disk phase—the origin of the
abnormally high levels of gas in its disk remains a mystery. The degree of axisymmetry of
the gas disk can provide clues as to the origin of the gas, which in turn allows new insight
into how gas-rich debris disks retain such large amounts of gas for so long. If the gas is ax-
isymmetric, it implies a steady-state evolution (consistent with primordial origin or ongoing
evaporation of cometary bodies), but if it is non-axisymmetric, it provides evidence for a recent
collision of Mars-sized planetary bodies or resonances induced by a giant planet. Because 49
Ceti is viewed at a high inclination relative to Earth, it is necessary to account for the view-
ing geometry by performing a “deprojection” of the interferometric data (or visibilities) that
assumes an underlying circular geometry before the axisymmetry can be assessed. The summer after my freshman year of college, I extended
the functionality of Professor Meredith Hughes’ old visibility deprojection code, rewriting it
in Python and updating it so that it can include spectral lines— previously it could only handle
continuum. We used this code to analyze Atacama Large Millimeter Array (ALMA) observa-
tions of the distribution of gas and dust in 49 Ceti’s disk, and found the gas to be predominantly
axisymmetric, lending credence to a steady-state evolution theory.

For more information on 49 Ceti, please take a look at my paper published in the Keck Northeast Astronomy Consortium's symposium proceedings, ["Searching for Non-Axisymmetry in 49 Ceti’s Unusual Gas-Rich Debris Disk."] (49Ceti_Paper.pdf)


###Deprojection Code
The deprojection code I wrote requires an input .uvf file containing disk data, and has several applications.

It can radially average the disk and compare it to a model, plotting flux density as a function of distance from the center of the disk. For example:
![Image](/49Ceti_continuum_dep_averaged.png)

It can also produce continuum images, channel maps, or moment maps. Moment map of 49 Ceti created by my code with a gaussian taper applied applied and plotted using code originally written by Jesse Lieman-Sifry:
![Image](/49Ceti_line_dep_taper_moment.png)
