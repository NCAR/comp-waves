
;PURPOSE: Subroutine to compute phase speed by cross correlating rows of the velocity
;         space-time diagram (which is input). The time series used for cross correlation 
;         is formed by collapsing all the rows of the diagram after shifting by an 
;         estimate of the phase speed. This estimate is obtained using only the central 
;         rows of the diagram.
;  
;
;
;
;  INPUTS: vel - filtered time-distance velocity map. Either pro or retrograde.
;          delx - pixel size in Mm
;          delt - cadence of observations
;
;  KEYWORDS:   debug - obvious...
;
;              ret   - set this if the computed speed is a retrograde
;                      speed (for a better initial guess for the fit).
;         
;              force_zero - force the fit to go through (0,0), i.e.
;                           through the reference pixel. By default, it
;                           it is turned off in the wave tracking code,
;                           because it gives higher errors even though
;                           the fit might still be good (i.e. same
;                           slope, but with an offset in x/y). 
;HISTORY: Created by ? ?
;         Corrected in Feb 2013 by C. Bethge. More accurate phase speeds and
;         errors now.
;         R Morton 11/2014 - Added code to deal with even length tracks.
;                            Added extra constraint for CC - max CC value has to be
;                            above probability of that of random signals.
;                            Commented code.
;                            Corrected bugs in phase speed error calculation
;                            Added new method for rejection of outliers
;                            Added correction to t-factor for phase speed error to account for varying track lengths
;         R Morton 05/2015   Fixed bugs in errors for phase speed fitting

FUNCTION compute_speed_newt, vel, delx, delt, debug=debug, ret=ret, force_zero=force_zero


s=size(vel)
nt=s(1)
npt=s(2)
IF npt mod 2 EQ 0 THEN npt=npt-1 ;deals with even length tracks
icent=npt/2


nlag=npt  ;number of points in cross correlation (make odd)
lag=indgen(nlag)-fix(nlag/2.)





;##############################################################
;Cross-correlate adjacent rows in time-distance map
;##############################################################

split=nt/25
cent=fltarr(npt,split)
spc=fltarr(split,2)
bnd=nt/split
spc[*,0]=bnd*findgen(split)
spc[*,1]=[spc[1:split-1]-1,nt-1]
atrous,vel,decomp=a,n_scales=1
vel=a[*,*,0] ;+smooth(a[*,*,1],3,/edge_truncate)

FOR i=-npt/2,npt/2 DO BEGIN
  FOR j=0,split-1 DO BEGIN
      in1=vel(spc[j,0]:spc[j,1],icent)
      in2=vel(spc[j,0]:spc[j,1],icent+i)
      cent(i+npt/2,j)=ccpeak(in1,in2,bnd/2)
      ;ccor=c_correlate(in1,in2,lag)
      ;m=max(ccor,imax)
      ;imax=1 > imax < (nlag-2)
      ;cent(i+n/2,j)=parabola(lag(imax-1:imax+1),1.-ccor(imax-1:imax+1))
  ENDFOR

ENDFOR

;##############################################################
;bootstrapping of means
;##############################################################
numbs=500
g=fix(randomu(systime,split,numbs)*(split))
centbs=fltarr(npt,split,numbs)
for i=0,split-1 do centbs[*,i,*]=cent[*,g[i,*]]

meanvar=(moment(mean(centbs,dim=2),dim=2) )[*,0:1]

meanvar[icent,1]=1.

;sometimes poor CC values which give 0 Stand Dev.
in=where(meanvar[*,1] eq 0)
IF n_elements(in) GT 0 THEN meanvar[in,*]=max(meanvar[*,1])


;##############################################################
;Fit line with constraints on fit
;The y intercept is constrained to go through 0 (.fixed)
;The gradient is confined to be between 0.001 and 20 (.limited & .limits)
;##############################################################
parinfo = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0.D]}, 2)
parinfo[1].limited[0] = 1 
parinfo[1].limited[1] = 1 
parinfo[0].fixed = 1

if keyword_set(ret) then begin
 parinfo[1].limits[0]  = -20D
 parinfo[1].limits[1]  = -0.001D
 parinfo[*].value      = [0D,-0.5D] 
endif else begin
 parinfo[1].limits[0]  = 0.001D
 parinfo[1].limits[1]  = 20D
 parinfo[*].value      = [0D,0.5D] 
endelse

;xtol limits accuracy of iterative procedure - here limited to 5dp - this is for speed
param = mpfitfun('lin_func_fit', lag, meanvar[*,0], sqrt(meanvar[*,1]), parinfo=parinfo, /quiet,bestnorm=bstn,perror=per,xtol=1d-5)
i_yfit = param[0] + lag*param[1]

IF keyword_set(debug) THEN BEGIN
   plot,lag,meanvar[*,0]
   oploterr,lag,meanvar[*,0],sqrt(meanvar[*,1])
   oplot,lag,i_yfit,linestyle=2
   ;pause

ENDIF


phase_speed = 1./param[1]*delx/delt 
speed_error=per[1]/param[1]^2*delx/delt 


return,[phase_speed,speed_error]


end
