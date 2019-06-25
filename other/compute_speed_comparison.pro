function compute_speed_comparison, vel, delx, delt, debug=debug, ret=ret, force_zero=force_zero
;  Subroutine to compute phase speed by cross correlating rows of the velocity
;  space-time diagram (which is input). The time series used for cross correlation is
;  formed by collapsing all the rows of the diagram after shifting by an estimate of the
;  phase speed. This estimate is obtained using only the central rows of the diagram.
;  Note, input delx in Mm!
;
;  corrected in Feb 2013 by C. Bethge. More accurate phase speeds and
;  errors now.
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


s=size(vel)
nt=s(1)
npt=s(2)
icent=npt/2

nlag=npt  ;number of points in cross correlation (make odd)
lag=indgen(nlag)-fix(nlag/2.)
n = 19 < npt  ;number of points along track for initial guess (make odd)
cent=fltarr(n)
for i=-n/2,n/2 do begin
  ccor=c_correlate(vel(*,icent),vel(*,icent+i),lag)
  m=max(ccor,imax)
  imax=1 > imax < (nlag-2)
  cent(i+n/2)=parabola(lag(imax-1:imax+1),1.-ccor(imax-1:imax+1))
;  cent(i+n/2)=total(lag(imax-1:imax+1)*ccor(imax-1:imax+1))/total(ccor(imax-1:imax+1))
endfor

c=[0.,median(deriv(cent))]   ;use robust estimate for mean shift;  collapse all time series along track
col=fltarr(nt)

for i=0,npt-1 do begin
  shift=c(1)*(i-npt/2)+c(0)
  col=col+interpolate(vel(*,i),findgen(nt)+shift)
endfor
col=col/float(npt)
;  use collapsed time series to perform cross correlation with all time series along track

cent=fltarr(npt)
dcent=fltarr(npt)+.1

if keyword_set(debug) then begin
 plot,lag,/nodata,xr=[-nlag/2.,nlag/2.],yr=[-npt/4.-1,npt/4.+1],ysty=1,xsty=1,xtit='Lag (Time Steps)',$
 	ytit='Correlation',chars=1.5
endif

;first iteration; for an initial guess to determine a useful
;range for searching the position of max correlation
for i=0,npt-1 do begin
  ccor=c_correlate(col,vel(*,i),lag)
  if keyword_set(debug) then begin
   oplot,lag,ccor-npt/4.+i*0.5, col=50
  endif
  m=max(ccor,imax)
  imax=1 > imax < (nlag-2)
  cent(i)=parabola(lag(imax-1:imax+1),1.-ccor(imax-1:imax+1))
endfor

shifts=deriv(cent)*delt/delx
init_guess=lag*median(shifts)
old_speeds=1./shifts
old_phase_speed=median(old_speeds)
old_errs=abs(old_speeds-old_phase_speed)
old_speed_error=median(old_errs)/sqrt(npt)

pdistance = fltarr(npt)
for ii=0,npt-1 do begin
 pdistance[ii] = pnt_line([lag[ii],cent[ii]], [lag[0],init_guess[0]], [lag[npt-1],init_guess[npt-1]])
endfor

good_p = where(pdistance le 2.*median(pdistance))
err = 1.

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

param = mpfitfun('lin_func_fit', lag[good_p], cent[good_p], err, parinfo=parinfo, /quiet)
i_yfit = param[0] + lag*param[1]

tstep_tolerance = 4
for i=0,npt-1 do begin
  ccor=c_correlate(col,vel(*,i),lag)
  lo_l = i_yfit[i]-tstep_tolerance+npt/2. ; lower limit
  up_l = i_yfit[i]+tstep_tolerance+npt/2. ; upper limit
  if (lo_l gt 1 and up_l lt nlag-2) then begin 
    m=max(ccor[lo_l:up_l],imax)
    imax=imax+i_yfit[i]-tstep_tolerance+npt/2.
    if keyword_set(debug) then plots, imax-nlag/2., m -npt/4. + i*0.5, psym = 4, color=100
    if keyword_set(debug) then plots, lo_l-nlag/2., ccor[lo_l] -npt/4. + i*0.5, psym = 5, color=250
    if keyword_set(debug) then plots, up_l-nlag/2., ccor[up_l] -npt/4. + i*0.5, psym = 5, color=250
  endif else begin
   m=max(ccor,imax)
   imax=1 > imax < (nlag-2)
   if keyword_set(debug) then plots, imax-nlag/2., m -npt/4. + i*0.5, psym = 4, color=185
  endelse
  cent(i)=parabola(lag(imax-1:imax+1),1.-ccor(imax-1:imax+1))
if keyword_set(debug) then plots, cent[i], m -npt/4. + i*0.5, psym = 5
endfor

data = cent*delt

xfit=(findgen(npt)-fix(npt/2))*delx
if keyword_set(debug) then begin
 plot,xfit,data,psym=1,xtit='Position Along Track (Mm)',ytit='Lag (s)',chars=1.5,$
     xr=[-npt/2,npt/2]*delx*1.5,xsty=1,yr=[-npt,npt]*4,ysty=1
endif

shifts=deriv(cent)*delt/delx
init_guess=xfit*median(shifts)

;first iteration
pdistance = fltarr(npt)
for ii=0,npt-1 do begin
 pdistance[ii] = pnt_line([xfit[ii],data[ii]], [xfit[0],init_guess[0]], [xfit[npt-1],init_guess[npt-1]])
endfor

good_p = where(pdistance le 2.*median(pdistance))
err = 1.

parinfo = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0.D]}, 2)
parinfo[1].limited[0] = 1 
parinfo[1].limited[1] = 1 
if keyword_set(force_zero) then begin
  parinfo[0].fixed = 1
endif

if keyword_set(ret) then begin
 parinfo[1].limits[0]  = -20D
 parinfo[1].limits[1]  = -0.001D
 parinfo[*].value      = [0D,-0.5D] 
endif else begin
 parinfo[1].limits[0]  = 0.001D
 parinfo[1].limits[1]  = 20D
 parinfo[*].value      = [0D,0.5D] 
endelse

param = mpfitfun('lin_func_fit', xfit[good_p], data[good_p], err, parinfo=parinfo, /quiet)
i_yfit = param[0] + xfit*param[1]

;second iteration
pdistance = fltarr(npt)
for ii=0,npt-1 do begin
 pdistance[ii] = pnt_line([xfit[ii],data[ii]], [xfit[0],i_yfit[0]], [xfit[npt-1],i_yfit[npt-1]])
endfor

good_p = where(pdistance le 10.*median(pdistance),comp=bad_p)
num_good_p = n_elements(good_p)
fin_val = where(finite(pdistance) eq 1,num_fin_val)
if num_fin_val eq n_elements(pdistance) then begin
 err = pdistance[good_p]
endif else begin
 err = 1.
endelse

param = mpfitfun('lin_func_fit', xfit[good_p], data[good_p], err, parinfo=parinfo, /quiet)
new_yfit = param[0] + xfit*param[1]

if num_good_p gt 0.3333*num_fin_val then begin
 phase_speed = 1./param[1] 
 x_m = (total(xfit[good_p])^2.)/float(num_good_p)
 s_xx = (total((xfit[good_p])^2.)-x_m)
; ~95% confidence level
; t_factor = 2.
; speed_error =((sqrt(total((data[good_p]-new_yfit[good_p])^2.)/(num_good_p-2)))/(sqrt(s_xx)))*2. 
; ~68% confidence level (different t-factors)
 t_factor = 1.02
 speed_error = ((sqrt(total((data[good_p]-new_yfit[good_p])^2.)/(float(num_good_p)-2.)))/(sqrt(s_xx)))*t_factor
endif else begin
 phase_speed = !Values.F_NAN
 speed_error = !Values.F_NAN
endelse

if keyword_set(debug) then begin
 if n_elements(bad_p) eq 1 then begin 
  if bad_p[0] ne -1 then begin
   oplot,xfit[bad_p],data[bad_p], color=230, psym=1
  endif
 endif else begin
   oplot,xfit[bad_p],data[bad_p], color=230, psym=1
 endelse
 oplot,xfit,new_yfit, color = 150
 oplot,xfit,init_guess, color = 185
endif

return,[phase_speed,speed_error,old_phase_speed,old_speed_error]
end
