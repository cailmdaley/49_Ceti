pro mfits,in=infile,in2=otherfile,ps=ps,vla=vla,block=block

; This procedure will read in a fits file (infile), and bin and radially 
; average the deprojected visibilities. 
; If you have a model you would like to compare it to, include it as otherfile.
; /ps will write out a postscript plot.
; use /vla if your data come from the VLA (slightly different uvfits 
;   convention). 
; /block will make interlocking bins
; NOTE: Please don't forget to include the proper position offsets in the
;   position section.  

if not keyword_set(infile) then begin
  print,'procedure mfits'
  print,'syntax: mfits,in=filename[,in2=filename,/ps,/vla,/block] '
  print,'please provide a filename.'
  return
endif

; adjustable parameters
nbins=17
; NOTE: PA is 180-(PA on sky)
pa=!dtor*double(35.0)			; NOTE: This PA came from Hamidouche et al 2006, NOT from Pietu et al who are just wrong (they quote 58 deg instead of 145)
incl=!dtor*double(38.0)			
; position offset in arcsec (see pos offset section below)
adrain=0.0
addecin=0.0


; reading in fits file
data=mrdfits(infile,0,header,/dscale,/unsigned,status=status)
freq=double(sxpar(header,'crval4'))
if keyword_set(uvf) then begin
freq=double(sxpar(header,'crval5'))
endif
if keyword_set(otherfile) then begin
  odata=mrdfits(otherfile,0,header2,/dscale,/unsigned,status=status)
  freq2=double(sxpar(header2,'crval4'))
  amp2=sqrt((odata.array[0])^2+(odata.array[1])^2)
endif

; rearranging resultant array of structures
if not keyword_set(vla) then begin
  amp=sqrt((data.array[0])^2+(data.array[1])^2)
  phase=atan(data.array[1]/data.array[0])
  real=data.array[0]
  imag=data.array[1]
  weight=data.array[2]
endif else begin
  q=data.array
  real=[reform(q(0,0,0,0,*)),reform(q(0,0,0,1,*))]
  imag=[reform(q(1,0,0,0,*)),reform(q(1,0,0,1,*))]
  weight=[reform(q(2,0,0,0,*)),reform(q(2,0,0,1,*))]
  amp=sqrt(real^2+imag^2)
  phase=atan(imag/real)
endelse
if keyword_set(otherfile) then begin
  amp2=sqrt((odata.array[0])^2+(odata.array[1])^2)
  phase2=atan(odata.array[1]/odata.array[0])
  real2=odata.array[0]
  imag2=odata.array[1]
  weight2=odata.array[2]
endif

; resolving quadrant ambiguities in phase
loc1=where((real lt 0d0) and (imag lt 0d0))
loc2=where((real lt 0d0) and (imag gt 0d0))
if loc1(0) ne -1 then phase(loc1)=phase(loc1)+!pi
if loc2(0) ne -1 then phase(loc2)=phase(loc2)+!pi

; calculating u-v distance of each point 
if not keyword_set(vla) then begin
  u=data.params[0]
  v=data.params[1]
endif else begin
  w=data.params
  u=[reform(w(0,*)),reform(w(0,*))]
  v=[reform(w(1,*)),reform(w(1,*))]
endelse
uvdist=sqrt(u^2+v^2)
urad=u*freq
vrad=v*freq

; finding deprojected r_uv
theta=atan(v/u)
tpa=theta-pa
da=uvdist*cos(tpa)
;*cos(incl)
dbd=uvdist*sin(tpa)/cos(incl)
ruv=sqrt(da^2+dbd^2)
if keyword_set(otherfile) then begin
  v2=odata.params[1]
  u2=odata.params[0]
  uvdist2=sqrt(u2^2+v2^2)
  theta2=atan(v2/u2)
  tpa2=theta2-pa
  da2=uvdist2*cos(tpa2)*cos(incl)
  dbd2=uvdist2*sin(tpa2)
  ruv2=sqrt(da2^2+dbd2^2)
endif

; ************** beginning of loop to fit for positions *****************
;  nposx=25		; must be odd
;  nposy=25		; ditto.
;  xinit=0.20
;  yinit=-0.23
;  xstep=0.01
;  ystep=0.01
;  flux=dblarr(nposx,nposy)
;  flux2=dblarr(nposx,nposy)
;  x=dindgen(nposx)*xstep+xinit-xstep*double(fix(nposx/2))
;  y=dindgen(nposy)*ystep+yinit-ystep*double(fix(nposy/2))
;  for k=0,nposx-1 do begin
;    for l=0,nposy-1 do begin
;      adra=x[k] & addec=y[l]


; calculating phase shifts due to pos offset
; Convention: this is the *shift* required; if the source center is -0.1,-0.1
; then you would choose adra=+0.1 & addec=+0.1
adra=adrain		; arcsec  MWC480 (SMA 345 GHz)
addec=addecin
dra=double(adra)*!pi/(3600d0*180d0)
ddec=double(addec)*!pi/(3600d0*180d0)
dphas=2d0*!pi*(urad*dra+vrad*ddec)
newphas=phase+dphas
newreal=amp*cos(newphas)
newimag=amp*sin(newphas)

if keyword_set(otherfile) then begin
; If there is a position offset in the comparison file, add it here
  adra2=0.0
  addec2=0.0
  dra2=double(adra2)*!pi/(3600d0*180d0)
  ddec2=double(addec2)*!pi/(3600d0*180d0)
  urad2=u2*freq2
  vrad2=v2*freq2
  dphas2=2d0*!pi*(urad2*dra2+vrad2*ddec2)
  phase2=phase2+dphas2
  real2=amp2*cos(phase2)
  imag2=amp2*sin(phase2)
endif

; setting up bins
uvmin=min(ruv) & uvmax=max(ruv)
uvrange=uvmax-uvmin
bmin=dindgen(nbins)*uvrange/double(nbins)+uvmin
bmax=(dindgen(nbins)+1d0)*uvrange/double(nbins)+uvmin
bmid=(bmin+bmax)/2d0
bimag=dblarr(nbins) & breal=dblarr(nbins) & binsize=dblarr(nbins)
bsd=dblarr(nbins) & bamp=dblarr(nbins) & brsd=dblarr(nbins) & bisd=dblarr(nbins)

; interlocking bins, if desired
if keyword_set(block) then begin
  nibins=nbins-1
  ibmin=bmid(0:nbins-2) & ibmax=bmid(1:nbins-1)
  ibmid=(ibmin+ibmax)/2d0
  ibimag=dblarr(nibins) & ibreal=dblarr(nibins) 
  ibsd=dblarr(nibins) & ibamp=dblarr(nibins) & ibrsd=dblarr(nibins)
  ibisd=dblarr(nibins)
endif

; binning and (vector) averaging
for i=0,nbins-1 do begin
  loc=where((ruv ge bmin[i]) and (ruv le bmax[i]))
  if (loc(0) eq -1) then begin
    bimag[i]=0d0 & breal[i]=0d0
    bsd[i]=0d0
  endif else begin
  ntot=double(n_elements(loc)-1)
  bimag[i]=total(newimag(loc)*weight(loc),/nan)/total(weight(loc))
  breal[i]=total(newreal(loc)*weight(loc),/nan)/total(weight(loc))
  bamp[i]=sqrt(bimag[i]^2+breal[i]^2)
  bsd[i]=sqrt(total(weight(loc)*(amp(loc)-bamp[i])^2)/total(weight(loc)))
  brsd[i]=sqrt(total(weight(loc)*(real(loc)-breal[i])^2)/(ntot*total(weight(loc))))
  bisd[i]=sqrt(total(weight(loc)*(imag(loc)-bimag[i])^2)/(ntot*total(weight(loc))))
  binsize[i]=double(n_elements(loc))
  endelse
  print,'number of elements in bin: ',binsize[i]
endfor

if keyword_set(block) then begin
for i=0,nibins-1 do begin
  loc=where((ruv ge ibmin[i]) and (ruv le ibmax[i]))
  if (loc(0) eq -1) then begin
    ibimag[i]=0d0 & ibreal[i]=0d0
    ibsd[i]=0d0
  endif else begin
    ntot=double(n_elements(loc)-1)
    ibimag[i]=total(newimag(loc)*weight(loc),/nan)/total(weight(loc))
    ibreal[i]=total(newreal(loc)*weight(loc),/nan)/total(weight(loc))
    ibamp[i]=sqrt(ibimag[i]^2+ibreal[i]^2)
    ibsd[i]=sqrt(total(weight(loc)*(amp(loc)-ibamp[i])^2)/total(weight(loc)))
    ibrsd[i]=sqrt(total(weight(loc)*(real(loc)-ibreal[i])^2)/(ntot*total(weight(loc))))
    ibisd[i]=sqrt(total(weight(loc)*(imag(loc)-ibimag[i])^2)/(ntot*total(weight(loc))))
  endelse
endfor
endif


; ***************** end of position loop ********************
;      flux[k,l]=mean(bamp)
;      flux[k,l]=total(abs(breal(7:10)))
;       flux[k,l]=total(brsd(5:10))
;       flux2[k,l]=total(bisd(5:10))
;      flux[k,l]=total(breal)
;      flux2[k,l]=total(abs(bimag))
;       flux[k,l]=total(abs(bimag(0:8)))
;       flux[k,l]=mean(abs(bimag))
;       flux2[k,l]=abs(stddev([bimag,ibimag]))
;      flux[k,l]=total(bisd(0:6))+total(brsd(0:6))+total(ibisd(0:6))+total(ibrsd(0:6))
;      flux[k,l]=abs(stddev([bimag(0:3),ibimag(0:2)]))
;      flux2[k,l]=abs(mean([bimag(0:3),ibimag(0:2)]))
; total(brsd(0:6))
;total(ibisd(0:6))+total(ibrsd(0:6))
;      flux[k,l]=breal(0)
;      flux[k,l]=bamp(0)
;    endfor
;  endfor
;  minmax=max(flux)
;  loc=where(flux eq minmax)
;  loc2=where(flux2 eq min(flux2))
;  yloc=fix(loc/nposx)
;  xloc=loc-nposx*yloc
;  yloc2=fix(loc2/nposx)
;  xloc2=loc2-nposx*yloc2
;  print,'max total flux is ',minmax
;  print,' and occurs at the position '
;  print,'delta  ra = ',x[xloc]
;  print,'delta dec = ',y[yloc]
;  print,'Criterion #2 min occurs at '
;  print,'delta ra  = ',x[xloc2]
;  print,'delta dec = ',y[yloc2]
;stop  

print,'total(bisd(1:4))=',total(bisd(1:4))
print,'total(brsd(1:4))=',total(brsd(1:4))


; plotting section
; ----------------
if keyword_set(ps) then begin
  set_plot,'ps',/interpolate
  device,filename='mfits.ps',bits_per_pixel=8 $
    ,xsize=10.0,ysize=7.5,xoff=0.5,yoff=10.5,/inches,/land,/color
endif


if keyword_set(otherfile) then begin
  ; setting up bins
  nbins2=15
  uvmin2=min(ruv2) & uvmax2=max(ruv2)
  uvrange2=uvmax2-uvmin2
  bmin2=dindgen(nbins2)*uvrange2/double(nbins2)+uvmin2
  bmax2=(dindgen(nbins2)+1d0)*uvrange2/double(nbins2)+uvmin2
  bmid2=(bmin2+bmax2)/2d0
  ; binning and (vector) averaging
  bimag2=dblarr(nbins2) & breal2=dblarr(nbins2)
  bsd2=dblarr(nbins2) & bamp2=dblarr(nbins2)
  for i=0,nbins2-1 do begin
    loc=where((ruv2 gt bmin2[i]) and (ruv2 lt bmax2[i]) and (newreal ne 0.0))
    if (loc(0) eq -1) then begin
      bimag2[i]=0d0 & breal2[i]=0d0
      bsd2[i]=0d0
    endif else begin
    bimag2[i]=total(imag2(loc)*weight2(loc))/total(weight2(loc))
    breal2[i]=total(real2(loc)*weight2(loc))/total(weight2(loc))
    bamp2[i]=sqrt(bimag2[i]^2+breal2[i]^2)
    bsd2[i]=sqrt(total(weight2(loc)*(amp2(loc)-bamp2[i])^2)/total(weight2))
    endelse
    endfor
endif 

; putting x-axis in klam
klam=bmid*freq*1.0e-3
if keyword_set(block) then iklam=ibmid*freq*1.0e-3
if keyword_set(otherfile) then klam2=bmid2*freq2*1.0e-3

  yratop=max(breal)-min(breal)
  yrabot=max(bimag)-min(bimag)
  xmax=max([klam])
!P.FONT=-1
q=findgen(17)*(!pi*2d0/16d0)
usersym,1.5*sin(q),1.5*cos(q),/fill
if keyword_set(xlim) then begin
  plot,klam,breal,xtitle='!12R!3!Duv!N (k!4k!3)',ytitle='real (Jy)', $
    /nodata,xrange=[bmid(0),xmax*1.1], $
    yrange=[min([breal,ibreal])-brsd(4),0.004],ystyle=1, $
    charthick=5,charsize=1.5,thick=3,xthick=6,ythick=6, $
    yminor=2,position=[0.20,0.8*yrabot/(yratop+yrabot)+0.15,0.95,0.95],xstyle=4
  axis,xaxis=1,xtickname=[" "," "," "," "," "," "],xthick=6
endif else begin
  plot,klam,breal,xtitle='!12R!3!Duv!N (k!4k!3)',ytitle='real (Jy)',/nodata, $
  charthick=5,charsize=2.0,thick=4,xthick=6,ythick=6, $
    yminor=2,position=[0.20,0.8*yrabot/(yratop+yrabot)+0.15,0.95,0.95],xstyle=5,$
    yrange=[min(breal)-yratop*0.1,max(breal)+yratop*0.1],xrange=[0.0,xmax*1.1],xminor=4,ystyle=1
  axis,xaxis=1,xtickname=[" "," "," "," "," "," "],xthick=6
endelse
oplot,klam,breal,psym=8,thick=6
merrplot,klam,breal-brsd,breal+brsd,width=0.01,errthick=6
usersym,1.5*sin(q),1.5*cos(q),thick=5
if keyword_set(block) then oplot,iklam,ibreal,psym=8  
if keyword_set(block) then merrplot,iklam,ibreal-ibrsd,ibreal+ibrsd,width=0.01,errthick=6
loadct,30,/silent
if keyword_set(otherfile) then oplot,klam2,breal2,linestyle=0,thick=4,color=100
loadct,0,/silent
oplot,[0,xmax*1.1],[0,0],linestyle=0,thick=4

plot,klam,bimag,xtitle='!12R!3!Duv!N (k!4k!3)',ytitle='imaginary (Jy)', $
  /nodata,xstyle=1,charthick=5,charsize=2.0,thick=3,xthick=6,ythick=6, $
  yminor=1,position=[0.20,0.15,0.95,0.8*yrabot/(yratop+yrabot)+0.15],/noerase, $
  yrange=[min(bimag)-yrabot*0.2,max(bimag)+yrabot*0.2],xminor=4,xrange=[0.0,xmax*1.1],ystyle=1
usersym,1.5*sin(q),1.5*cos(q),/fill
oplot,klam,bimag,psym=8,thick=6
merrplot,klam,bimag-bisd,bimag+bisd,width=0.01,errthick=6
usersym,1.5*sin(q),1.5*cos(q),thick=5
if keyword_set(block) then oplot,iklam,ibimag,psym=8 
if keyword_set(block) then merrplot,iklam,ibimag+ibisd,ibimag-ibrsd,width=0.01,errthick=6
oplot,[0,xmax*1.1],[0,0],linestyle=0,thick=4
radist=51.0*3600.0/(klam*1000.0*!dtor)

if keyword_set(ps) then begin
  device,/close
  set_plot,'x'
  !p.font=-1
  print,"created postscript file: mfits.ps"
endif


return
end


; $Id: errplot.pro,v 1.16 2004/01/21 15:54:52 scottm Exp $
;
; Copyright (c) 1983-2004, Research Systems, Inc.  All rights reserved.
;   Unauthorized reproduction prohibited.
;

Pro merrplot, X, Low, High, Width = width, $
    DEVICE=device, $   ; swallow this keyword so user can't set it
    NOCLIP=noclipIn, $ ; we use a different default than PLOTS
    _REF_EXTRA=_extra,errthick=errthick
;+
; NAME:
;   ERRPLOT
;
; PURPOSE:
;   Plot error bars over a previously drawn plot.
;
; CATEGORY:
;   J6 - plotting, graphics, one dimensional.
;
; CALLING SEQUENCE:
;   ERRPLOT, Low, High	;X axis = point number.
;
;   ERRPLOT, X, Low, High   ;To explicitly specify abscissae.
;
; INPUTS:
;   Low:    A vector of lower estimates, equal to data - error.
;   High:   A vector of upper estimates, equal to data + error.
;
; OPTIONAL INPUT PARAMETERS:
;   X:	A vector containing the abscissae.
;
; KEYWORD Parameters:
;   WIDTH:  The width of the error bars, in units of the width of
;   the plot area.  The default is 1% of plot width.
;
;   All keywords to PLOTS are also accepted.
;
; OUTPUTS:
;   None.
;
; COMMON BLOCKS:
;   None.
;
; SIDE EFFECTS:
;   An overplot is produced.
;
; RESTRICTIONS:
;   Logarithmic restriction removed.
;
; PROCEDURE:
;   Error bars are drawn for each element.
;
; EXAMPLES:
;   To plot symmetrical error bars where Y = data values and
;   ERR = symmetrical error estimates, enter:
;
;	PLOT, Y		;Plot data
;	ERRPLOT, Y-ERR, Y+ERR	;Overplot error bars.
;
;   If error estimates are non-symetrical, enter:
;
;	PLOT,Y
;	ERRPLOT, Upper, Lower	;Where Upper & Lower are bounds.
;
;   To plot versus a vector of abscissae:
;
;	PLOT, X, Y	  ;Plot data (X versus Y).
;	ERRPLOT, X, Y-ERR, Y+ERR  ;Overplot error estimates.
;
; MODIFICATION HISTORY:
;   DMS, RSI, June, 1983.
;
;   Joe Zawodney, LASP, Univ of Colo., March, 1986. Removed logarithmic
;   restriction.
;
;   DMS, March, 1989.  Modified for Unix IDL.
;	KDB, March, 1997.  Modified to used !p.noclip
;	RJF, Nov, 1997.	   Removed unnecessary print statement
;		Disable and re-enable the symbols for the bars
;   DMS, Dec, 1998.    Use device coordinates.	Cleaned up logic.
;   CT, RSI, Jan 2001: Add _REF_EXTRA to pass keywords to PLOTS.
;-
on_error,2			;Return to caller if an error occurs
if n_params(0) eq 3 then begin	;X specified?
    up = high
    down = low
    xx = x
endif else begin		;Only 2 params
    up = x
    down = low
    xx=findgen(n_elements(up))	;make our own x
endelse

w = ((n_elements(width) eq 0) ? 0.01 : width) * $ ;Width of error bars
  (!x.window[1] - !x.window[0]) * !d.x_size * 0.5
n = n_elements(up) < n_elements(down) < n_elements(xx) ;# of pnts

; If user hasn't set NOCLIP, follow what is in !P.
; This is different than PLOTS, whose default is always NOCLIP=1.
noclip = (N_ELEMENTS(noclipIn) gt 0) ? noclipIn : !P.NOCLIP

for i=0,n-1 do begin		;do each point.
    xy0 = convert_coord(xx[i], down[i], /DATA, /TO_DEVICE) ;get device coords
    xy1 = convert_coord(xx[i], up[i], /DATA, /TO_DEVICE)
    plots, [xy0[0] + [-w, w,0], xy1[0] + [0, -w, w]], $
	[replicate(xy0[1],3), replicate(xy1[1],3)], $
	/DEVICE, NOCLIP=noclip, _STRICT_EXTRA=_extra,thick=errthick
endfor
end


