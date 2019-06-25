pro wave_tracking, date, nwl=nwl,                           $
                         own_data=own_data,                 $
                         compute_speeds=compute_speeds,     $
                         cross_corr=cross_corr,             $
                         plot_maps=plot_maps,               $
                         save_coh = save_coh,               $
                         name_addon = name_addon,           $
                         wr_speeds = wr_speeds,             $
                         choose_corr_box = choose_corr_box, $
                         debug = debug, altangle=altangle,  $
                         freq0=freq0

;  Procedure to determine the wave propagation angle from the cross
;  correlation of the velocity time series of a test point with
;  adjacent time series. The computation is performed in the Fourier
;  domain. Only contiguous points having coherence greater than a
;  threshold (typically 0.5) are used. SSWIDL is needed.
;
;  The default is to use CoMP L2 data, in which case you just need to
;  provide a date (and the correct inpath/outpath); the code (ideally)
;  does the rest. You can also provide your own data (see keyword
;  "own_data"). The code starts looking for data early in the day,
;  which is usually the best due to the better seeing conditions. The
;  code then looks for the longest contiguous time sequence, and
;  interpolates over data gaps of one missing file, which is quite
;  often the case. The data is filtered, the default is a central
;  filter position of 3.5 mHz and filter width of 1.5 mHz. Change the
;  variables freq0 and fwidth if you want to change the filtering. 
;  Also, check if the default values of lower_r and upper_r (lower and
;  upper radius of the field-of-view for masking) are suitable for
;  your data.
;  
;  INPUT:      date            - STRING in the format yyyymmdd. 
;                                Will work on HAO's machines or if you
;                                provide CoMP L2 data in your own folder
;                                (change inpath in this case). 
;  
;  KEYWORDS:   own_data        - STRING.
;                                You can provide your own Doppler velocities
;                                here. Set own_data to a sav-file containing
;                                the Doppler cube, i.e. 
;                                own_data='/home/me/doppler.sav'. In this
;                                case, make sure to rename the Doppler cube
;                                to cube_v and to set a couple of variables
;                                (norm_cadence, arcsec_per_pixel) further
;                                down in the code for a correct computation. 
;                                You will still have to provide date as input.
;
;              nwl             - INTEGER, 3 or 5. 
;                                Number of wavelengths (3pt or 5 pt data).
;                                Default is 3. 
;
;              cross_corr      - Remove (some) jitter from the CoMP data. If
;                                you do this, make sure to set the debug
;                                keyword for the first run, which will show
;                                you the box for which the cross-correlation 
;                                is computed. If the default box does not look 
;                                correct, set the values x1, x2, y1, y2
;                                manually below.  
;
;              choose_corr_box - Select the cross-correlation box with
;                                the mouse.
;
;              compute_speeds  - compute phase speeds after the
;                                computation of the wave propagation
;                                angle. 
;
;              plot_maps       - plot ps-files with all quantities after the
;                                computation.
;                              
;              name_addon      - STRING.
;                                Set this to a string that is added to the
;                                name of the sav-files and plots, to
;                                distinguish between different runs.                
;
;              save_coh        - Create a sav-file containing a measure of
;                                the coherence for each pixel. This is done
;                                by fitting an ellipse to the coherence
;                                island. The coherence for a pixel is good 
;                                when the ellipse is long and slim, i.e.
;                                the normalized ratio of the major and minor
;                                axis of the ellipse is a measure of the
;                                coherence. The sav-file will contain an
;                                array coh_measure, which is just the major
;                                (coh_measure[*,*,0]) and minor
;                                (coh_measure[*,*,1]) axis of the ellipse. 
;                                One possible way to compute a quantity for 
;                                the coherence quality would then be
;  coherence_quality = (coh_measure[*,*,0]-coh_measure[*,*,1])/(coh_measure[*,*,0]+coh_measure[*,*,1])
;
;             debug            - set this if you want to visually
;                                check what's going on. *MUCH MUCH*
;                                slower computation though, so don't
;                                use this as the default.
;                           
;             wr_speeds        - This is mainly intended for checking the
;                                computed speeds afterwards. It will write
;                                out *very* large arrays, so don't
;                                do this on a computer with insufficient
;                                memory or disk space, especially in
;                                combination with the save_coh keyword. Not
;                                really needed for everyday use, mostly for
;                                debugging. 
;
;              altangle - Use alternative estimate for wave angle based on 
;                         sixlin & choosing 


; use 3 pt data by default
if n_elements(nwl) eq 0 then nwl = 3
nwlst = strcompress(string(nwl),/rem)
ans   = ' '
torad = !pi/180.
if keyword_set(wr_speeds) then wr_speeds = 1 else wr_speeds = 0

;================================================================
;*** Here the details of the wave tracking need to be defined ***
;================================================================

; paths where the CoMP data is located and
; where the .sav-files should be written to

;IF keyword_set(cross_corr) THEN &
inpath  = 'data/CoMP/'+strmid(date,8,4,/reverse)+'/'+strmid(date,3,2,/reverse)+'/'+strmid(date,1,2,/reverse)+'/'
outpath = 'data/CoMP/wave_tracking_output/'+date+'/'


if file_test(outpath,/DIR) eq 0 then spawn, 'mkdir '+outpath
if keyword_set(plot_maps) then begin
 if file_test(outpath+'plots/',/DIR) eq 0 then spawn, 'mkdir '+outpath+'plots/'
endif

; values for the finding the wave propagation angles
dx      = 20      ;
nbox    = 2*dx+1  ;box size in pixels
nsmooth = 15      ;smoothing width for cross spectra (nominally 11)
limit   = 0.5     ;set coherence limit for acceptance
npix    = 10      ;minimum number of coherent pixels

IF n_elements(freq0) LT 1 THEN BEGIN
   ; central frequency and filter width
   ; for the frequency filtering in Hz
   ;freq0  = 0.0015
   freq0  = 0.0035
   ;freq0  = 0.0055
ENDIF

fwidth = 0.0015

;solar radius in Mm
rsun_mm = 696.342

;addon to the names of the .sav-files
if n_elements(name_addon) eq 0 then begin
 name_addon = ''
endif else begin
 name_addon = name_addon+'_'
endelse

case freq0 of
 0.0015: freqst = '1.5mHz'
 0.0035: freqst = '3.5mHz'
 0.0055: freqst = '5.5mHz'
 else: begin
         freqst = 'manual_freq_setting'
       end
endcase

; *** Default values for a date that has not ***
; *** been checked before and is not defined ***
; *** below in the case statement            ***

; values for cross-correlation
; x1,x2,y1,y2 - coordinates of the
; box for the cross-correlation
x1 = 45
x2 = 75
y1 = 230
y2 = 390

; For starting with a later file - or not (start_file = 0).
; Useful if the first (couple of) file(s) are not very good.
start_file = 0

; Number of files to process. 0 means: let the 
; code choose the highest number possible. Useful 
; for manually cutting out bad data at the end.
num_files_proc = 0

; Upper and lower limit for the radius
; !be aware that the maximum offset of the cross-correlation
; (if performed) is added to this value!
lower_r = 228.
upper_r = 288.

; threshold for the cross-correlation precision
ccthresh = 0.002

; maximum iteration steps for the cross-correlation
maxiter = 200. 

;==========================================
;manual processing values for certain dates
;==========================================
case date of
 ;example for a manual setting for a specific date
 '20130213': begin
              x1 = 545
              x2 = 580
              y1 = 375 ;236
              y2 = 410 ;270
              start_file = 0
              num_files_proc = 178
              lower_r = 228.
              upper_r = 288.
              ;maxiter = 300.
              ccthresh = 0.2
             end

  '20120327': begin
              x1 = 55
	      x2 = 85
	      y1 = 220 ;236
	      y2 = 260 ;270
	      start_file = 0
	      num_files_proc = 164
	      lower_r = 228.
	      upper_r = 288.
	      ;maxiter = 300.
	      ccthresh = 0.04
              end

  '20120410': begin
              x1 = 535
	      x2 = 565
	      y1 = 390
	      y2 = 430 
	      start_file = 0
	      num_files_proc = 300
	      lower_r = 228.
	      upper_r = 288.
	      ;maxiter = 300.
	      ccthresh = 0.04
              end

  else:      begin
               message, /cont, 'Manual settings for this date were not found, using default values.'
             end
endcase
box = [x1,x2,y1,y2]

if keyword_set(debug) then ccdebug = 1


;====== if keyword own_data is set  =================================
; skip the preparation of the data if someone provides his own
; velocity data cube (cube_v) - in this case, a few additional 
; variables need to be defined here: 
; norm_cadence = the cadence of cube_v in seconds
; arcsec_per_pixel = usually 4.46 for CoMP data
if n_elements(own_data) ne 0 then begin
 restore, own_data
 ;cube_v = velocity
 s=size(cube_v)
 nx=s[1] 
 ny=s[2]
 nt=s[3]
 
 norm_cadence     = 30.
 arcsec_per_pixel = 4.46
 pix              = arcsec_per_pixel*0.725   
                    ;pixels in Mm - assume 725 km/arcsec (1 AU distance). 
                    ;correct this value if necessary -
                    ;important for correct phase speeds!!
 goto, external_data
endif
;====== end if keyword own_data is set =============================

;================================================================
;***      End of defining the wave tracking parameters        ***
;================================================================


flist   = file_list(inpath,'*.comp.1074.dynamics.'+nwlst+'.fts.gz')
if start_file ne 0 then begin
 flist = flist[start_file:*]
endif
fnumber = n_elements(flist)
cadence = intarr(fnumber-1)
hdr     = headfits(flist[0],ext=1)
nx      = sxpar(hdr,'NAXIS1')
ny      = sxpar(hdr,'NAXIS2')

;================================================================
;compute cadence
;================================================================
for ii=0,fnumber-2 do begin
 hdr       = headfits(flist[ii])
 date_obs  = sxpar(hdr,'DATE-OBS')
 time_obs  = sxpar(hdr,'TIME-OBS')
 new_time  = anytim2tai(date_obs+' '+time_obs) 
 if ii gt 0 then begin
  cadence[ii-1] = fix(new_time-old_time)
  
 endif
 old_time = new_time
  
endfor

IF keyword_set(debug) THEN plot,cadence,xtitle='Frame No.',ytitle='Cadence (s)'
norm_cadence = median(cadence)

;================================================================
; look for the start of the normal cadence in case the sequence 
; starts with files with an abnormal cadence
;================================================================
seq_start = 0
seq_end   = fnumber-1
for ii=0,fnumber-2 do begin
 if cadence[ii] ne norm_cadence then begin
   goto, look_for_start
 endif else begin
   seq_start = ii
   goto, found_start
 endelse
look_for_start:
endfor
found_start:

;================================================================
; interpolate over gaps of one missing file
; (happens quite often), otherwise don't bother
;================================================================
itrplt_counter = 0.
itrplt_pos = intarr(1)
for ii=seq_start,fnumber-2 do begin
    if cadence[ii] ne norm_cadence then begin
         if cadence[ii] eq 2.*norm_cadence then begin
              if ii lt fnumber-3 then begin
                    if cadence[ii+1] eq norm_cadence then begin 
                       ; save interpolation positions   
                       itrplt_counter = itrplt_counter+1
                       itrplt_pos = intarr(itrplt_counter)
                       if itrplt_counter gt 1 then begin
                          itrplt_pos[0:itrplt_counter-2] = old_itrplt_pos
                          itrplt_pos[itrplt_counter-1] = ii
                       endif else begin
                          itrplt_pos[0] = ii
                       endelse
                       old_itrplt_pos = itrplt_pos
                    endif else begin
                       ; don't bother
                       seq_end = ii
                       goto, found_end
                    endelse 
               endif 
         endif else begin
            seq_end = ii 
            goto, found_end
        endelse 
    endif
endfor
found_end:



;================================================================
; cut out contiguous file list and prepare
; the interpolated data cube
;================================================================
print,'Start sequence',seq_start,' End sequence',seq_end
flist = flist[seq_start:seq_end]
itrplt_pos = itrplt_pos-seq_start
nt = seq_end-seq_start+1+itrplt_counter
cube_i = fltarr(nx,ny,nt)
cube_v = fltarr(nx,ny,nt)
cube_w = fltarr(nx,ny,nt)

pix = (rsun_mm/sxpar(headfits(flist[0]),'RSUN'))*sxpar(headfits(flist[0]),'CDELT1') ;Mm
;time_test = strarr(nt)+'ZERO' ; (in)sanity check

;================================================================
; fill cube_v with data and interpolate if necessary
;================================================================
if itrplt_counter gt 0 then begin
   if itrplt_counter gt 1 then begin
      message, /cont, strcompress(string(fix(itrplt_counter)),/rem)+' data gaps interpolated.'
   endif else begin
      message, /cont, strcompress(string(fix(itrplt_counter)),/rem)+' data gap interpolated.'
   endelse
   start_idx = 0
   for ii=0,itrplt_counter-1 do begin
       for jj=start_idx,itrplt_pos[ii]+ii do begin
           cube_i[*,*,jj] = readfits(flist[jj-ii],ext=1,/silent)
           cube_v[*,*,jj] = readfits(flist[jj-ii],ext=3,/silent)
           cube_w[*,*,jj] = readfits(flist[jj-ii],ext=4,/silent)
           ;   time_test[jj]  = sxpar(headfits(flist[jj-ii]),'TIME-OBS')
       endfor
       ;  print, itrplt_pos[ii]
       ;  time_test[itrplt_pos[ii]+ii+1]  = 'In between '+sxpar(headfits(flist[itrplt_pos[ii]]),'TIME-OBS')+$
       ;                                    ' and '+sxpar(headfits(flist[itrplt_pos[ii]+1]),'TIME-OBS')
       cube_i[*,*,itrplt_pos[ii]+ii+1] = (readfits(flist[itrplt_pos[ii]],ext=1,/silent)+$
                                     readfits(flist[itrplt_pos[ii]+1],ext=1,/silent))/2.
       cube_v[*,*,itrplt_pos[ii]+ii+1] = (readfits(flist[itrplt_pos[ii]],ext=3,/silent)+$
                                     readfits(flist[itrplt_pos[ii]+1],ext=3,/silent))/2.
       cube_w[*,*,itrplt_pos[ii]+ii+1] = (readfits(flist[itrplt_pos[ii]],ext=4,/silent)+$
                                     readfits(flist[itrplt_pos[ii]+1],ext=4,/silent))/2.
       start_idx = itrplt_pos[ii]+ii+2
   endfor
   for jj=start_idx,nt-1 do begin
       cube_i[*,*,jj] = readfits(flist[jj-itrplt_counter],ext=1,/silent)
       cube_v[*,*,jj] = readfits(flist[jj-itrplt_counter],ext=3,/silent)
       cube_w[*,*,jj] = readfits(flist[jj-itrplt_counter],ext=4,/silent)
       ;  time_test[jj]  = sxpar(headfits(flist[jj-itrplt_counter]),'TIME-OBS')
   endfor
endif else begin
 ; in case no interpolation is necessary
 for ii=0,nt-1 do begin
   ;  print, sxpar(headfits(flist[ii]),'TIME-OBS')
   cube_i[*,*,ii] = readfits(flist[ii],ext=1,/silent)
   cube_v[*,*,ii] = readfits(flist[ii],ext=3,/silent)
   cube_w[*,*,ii] = readfits(flist[ii],ext=4,/silent)
 endfor
endelse

; manually exclude bad data
if num_files_proc ne 0 then begin
  cube_i = reform(cube_i[*,*,0:num_files_proc])
  cube_v = reform(cube_v[*,*,0:num_files_proc])
  cube_w = reform(cube_w[*,*,0:num_files_proc])
endif

external_data:

;================================================================
; cross-correlate the data and align it
;================================================================
if keyword_set(cross_corr) then begin
  if keyword_set(choose_corr_box) then begin
     ONCE_MORE:
     x2 = -1
     y2 = -1
     window, /free, xs=620, ys=620
     win_nr = !d.window
     aia_lct, wave=193, /load
     tvscl, reform(cube_i[*,*,0]) > 0 < 25
     message, /cont, 'Choose first corner'
     cursor, x1, y1, /down, /device    
     message, /cont, 'Choose second corner'

     IF !MOUSE.BUTTON EQ 1 THEN BEGIN 
        CURSOR, do_nothing_x, do_nothing_y, /UP, /WAIT
        !MOUSE.BUTTON = 0  
     ENDIF
     while !Mouse.Button NE 1 do begin
        cursor, tmp_x2, tmp_y2, /change, /device
        tvscl, reform(cube_i[*,*,0]) > 0 < 25
        plots, [x1,tmp_x2], [y1,y1], /device, color=255 
        plots, [x1,tmp_x2], [tmp_y2,tmp_y2], /device, color=255  
        plots, [x1,x1], [y1,tmp_y2], /device, color=255 
        plots, [tmp_x2,tmp_x2], [y1,tmp_y2], /device, color=255 
     endwhile

     wdelete, win_nr
     if x1 EQ tmp_x2 OR y2 EQ tmp_y2 then begin
       message, /cont, 'Box must be at least 1x1 pixel' 
       goto, once_more
     endif
     if x1 GT tmp_x2 then begin
       x2 = x1
       x1 = tmp_x2
     endif else begin
       x2 = tmp_x2
     endelse
     if y1 GT tmp_y2 then begin
       y2 = y1
       y1 = tmp_y2
     endif else begin
       y2 = tmp_y2
     endelse
     message, /cont, 'Chosen box: x1='+strcompress(string(x1),/rem)+', x2='+strcompress(string(x2),/rem)+', y1='+$
              strcompress(string(y1),/rem)+', y2='+strcompress(string(y2),/rem)
     box = [x1,x2,y1,y2]
  endif
  
  dejitter,cube_i,cube_v,cube_w,box=box,ccthresh=ccthresh,maxiter=maxiter,ccdebug=ccdebug,maxoffset=maxoffset
  
endif
save, cube_i, cube_v, cube_w, filename=outpath+'cube_ivw_'+date+'.sav'
;stop
if n_elements(maxoffset) eq 0 then maxoffset = 0

nt = (size(cube_v))[3]
message, /cont, 'Total # of files: '+strcompress(string(nt,format='(F5.0)'),/rem)

;================================================================
;FFT of velocity data cube
;================================================================
nspec = nt/2
spec=complexarr(nx,ny,nspec,/nozero)
for iy=0,ny-1 do for ix=0,nx-1 do begin
  d=cube_v[ix,iy,*] ;-smooth(cube_v[ix,iy,*],40,/edge_truncate)
  sp=fft(d)
  spec[ix,iy,*]=sp[0:nspec-1]
endfor

s=size(spec)

freq   = findgen(nspec)/(float(nspec*2)*norm_cadence)
nfreq  = n_elements(freq)
filter = exp( -(freq-freq0)^2/fwidth^2 )
filter(0) = 0.                      ; set dc to zero
filter(where (freq lt .001)) = 0.   ; set low frequencies to zero
filter=filter/total(filter)

loadct, 39, /silent

; create mask
; mask out pixels to invert
mask=intarr(nx,ny)
x=rebin(findgen(ny)-(nx/2.)+0.5,nx,ny)
y=transpose(x)
r=sqrt(x^2+y^2)
good=where(r ge lower_r+maxoffset and r le upper_r)
mask(good)=1
if not keyword_set(debug) then begin
 pixcounter = 0L
 old_perc_proc = 0
 tot_pix = long(n_elements(good))
endif

; set x and y limits for map
; entire map
xstart=0
xend=nx-1
ystart=0
yend=ny-1

; open windows
f=8
if keyword_set(debug) then begin
 window,1,xs=nbox*f,ys=nbox*f,xpos=0,ypos=0
 window,4,xs=nx,ys=ny,retain=2,xpos=400+nx+20,ypos=0, title=date
 window,5,xs=nx,ys=ny,retain=2,xpos=400,ypos=ny+50
 window,6,xs=512,ys=512,xpos=400,ypos=800
 device,decomposed=0
endif

xbox=rebin(findgen(nbox)-dx,nbox,nbox)  ;coordinates of points in box
ybox=transpose(xbox)

wave_angle=fltarr(nx,ny)
angle_error=fltarr(nx,ny)
if keyword_set(save_coh) then coh_measure = fltarr(nx,ny,2)

;**********!!*!*!*!**!*!*!*!*!*!*!*!*!*!
;**********!!*!*!*!**!*!*!*!*!*!*!*!*!*!
;**********!!*!*!*!**!*!*!*!*!*!*!*!*!*!
;**********!!*!*!*!**!*!*!*!*!*!*!*!*!*!
;openw,lun,'data/comp/test.txt',/get_lun
angle_error2=fltarr(nx,ny)
coh_measure=fltarr(nx,ny,2)

message, /cont,'Starting processing...'
for iy=ystart,yend do for ix=xstart,xend do if mask(ix,iy) eq 1 then begin

  if not keyword_set(debug) then begin
   perc_proc = round((float(pixcounter)/float(tot_pix))*100.)
   pixcounter = pixcounter+1
   if perc_proc gt old_perc_proc then begin
    if strmid(getenv('TERM'),0,5) eq 'xterm' then begin
     progress_comp, perc_proc, 100, msg='Computing wave propagation angles'
    endif else begin
     message, /cont, strcompress(string(perc_proc,format='(I3)'),/rem)+'% processed.'
    endelse
    old_perc_proc = perc_proc
   endif
  endif

  ;================================================================
  ;  compute coherence for contiguous pixels with coherence greater than limit
  ;================================================================

  coh=fltarr(nbox,nbox)           ;initialize coherence array to zero
  coh_mask=fltarr(nbox,nbox)      ;initialize mask of good points in coherence

  spec1=spec[ix,iy,*]
  g1=real_part(spec1*conj(spec1))
  g1=smooth(g1,nsmooth)

  coh(nbox/2,nbox/2)=1.
  coh_mask(nbox/2,nbox/2)=1.

  count=1
  ncoh=1
  while count gt 0 and ncoh lt 100 do begin
    cont=where(coh_mask eq 1.)
    cnew=[cont-1,cont+1,cont+nbox,cont-nbox]  ;test pixels adjacent to contiguous pixels

    new_mask=fltarr(nbox,nbox)   ;identify new contiguous pixels
    new_mask(cnew)=1.
    new_mask=new_mask*(1.-coh_mask)

    new=where(new_mask gt 0.,new_count)

    count=0
    for ic=0,new_count-1 do begin
      icx=new(ic) mod nbox
      icy=fix(new(ic)/nbox)

      ii=ix-nbox/2+icx
      jj=iy-nbox/2+icy

      if ii lt 0 or ii gt nx-1 or jj lt 0 or jj gt ny-1 then coh(icx,icy)=0. else begin
        if coh(icx,icy) eq 0. then begin
          spec2=spec(ii,jj,*)
          cspec=spec1*conj(spec2)
          cspec=smooth(cspec,nsmooth)
          g2=real_part(spec2*conj(spec2))
          g2=smooth(g2,nsmooth)

          coh(icx,icy)=mask(ii,jj)*total( filter*( abs(cspec)/sqrt(g1*g2) ) ) 
          ;see McIntosh et al. 2008. Sol. Phys.
        endif

        if coh(icx,icy) gt limit then begin
          coh_mask(icx,icy)=1.
          count=count+1
        endif

      endelse

    endfor

   good=where(coh_mask eq 1.,ncoh)

   endwhile

   if keyword_set(debug) then begin
    wset,1
    tv,bytscl(rebin(coh*coh_mask,f*nbox,f*nbox,/sample),0.,1.)
   endif   

   ;if (keyword_set(save_coh) and (ncoh gt npix)) then begin
     island = coh*coh_mask
     ellipsepts = fit_ellipse(where(island ne 0.),axes=axes,xsize=nbox,ysize=nbox)
     coh_measure[ix,iy,*] = axes
     if ((size(ellipsepts))[0] gt 0) and keyword_set(debug) then plots, f*ellipsepts, color=254, /device
   ;endif

;================================================================
;  find angle of coherence island if enough good points
;================================================================

  if ncoh gt npix then begin
   if keyword_set(debug) then begin
    wset,6
    plot,xbox(good),ybox(good),psym=6,xr=[-dx,dx],yr=[-dx,dx]
   endif
 
    weight=coh(good)^2
    theta=min_dist_fit(xbox(good),ybox(good),error=error) ; analytical solution to the minimum distance problem
    ttest=min_test(xbox(good),ybox(good),error=err2)
 
    torad=!pi/180.
   ; printf,lun,theta,error
   ; printf,lun,'n',(ttest)/(torad),atan(err2)/(torad)
    ;theta= (ttest)/(torad)
    angle_error2(ix,iy)=atan(err2)/(torad)

    IF keyword_set(altangle) THEN BEGIN

        xe=xbox(good) & ye=ybox(good)
        sxy=total((xe-mean(xe))*(ye-mean(ye)))

        if round(sxy) ne 0. then begin
    
           sixlin,xbox(good),ybox(good),a,siga,b,sigb
    
           ;printf,lun,atan(theta*(!pi/180.)),atan(error*(!pi/180.))
           ;printf,lun,(b),(sigb)
           torad=!pi/180.
           theta=[atan(b[2])/(torad)]
           error=[atan(sigb[2])/(torad)]
         endif
      ENDIF
   
    if keyword_set(debug) then print,theta,error

    wave_angle(ix,iy)=theta
    angle_error(ix,iy)=error

    xfit=findgen(2*dx+1)-dx
    coef=[0.,tan(theta*torad)]
    yfit=poly(xfit,coef)  ;C0 + C1*x + C2*x^2
 
   if keyword_set(debug) then begin
    oplot,xfit,yfit
    wset,4
    tv,bytscl(wave_angle,-90,90)
    wset,5
    tv,bytscl(angle_error,0,10)
    wait,0.005
   endif
  endif
  
end

;**********!!*!*!*!**!*!*!*!*!*!*!*!*!*!
;**********!!*!*!*!**!*!*!*!*!*!*!*!*!*!
;**********!!*!*!*!**!*!*!*!*!*!*!*!*!*!
;free_lun,lun

;set_plot,'ps'
;device,/color,/encapsul,filename=outpath+'err_test.eps'
;!p.multi=[0,2,1]

;in=where(angle_error2 gt 0)
;in2=where(angle_error gt 0)
;ax1=coh_measure[*,*,0]
;ax2=coh_measure[*,*,1]

;plot,angle_error2[in],ax1[in]/ax2[in],psym=1,xtitle='angle error',ytitle='Ratio coherence length'
;plot,angle_error[in2],ax1[in2]/ax2[in2],psym=1,xtitle='angle error',ytitle='Ratio coherence length'

;device,/close
set_plot,'x'
;!p.multi=0


if keyword_set(save_coh) then save, coh_measure, file=outpath+'coherence_measure_'+name_addon+date+'_'+nwlst+$
                                    '_'+freqst+'.sav'

xscale = pix
save,wave_angle,angle_error,xscale,norm_cadence,mask,file=outpath+'wave_angle_'+name_addon+date+'_'+nwlst+'_'+freqst+'.sav'

if keyword_set(compute_speeds) then begin
 space_time_run, date, cube_v, wave_angle, mask, xscale, lower_r+maxoffset, $
                 norm_cadence, outpath, nwlst, freqst, name_addon, debug=debug, wr_speeds=wr_speeds
endif
if keyword_set(plot_maps) then begin
 if keyword_set(save_coh) then begin
  plot_wave_maps, date, nwl, freqst, outpath, name_addon, /cohm
 endif else begin
  plot_wave_maps, date, nwl, freqst, outpath, name_addon
 endelse
endif

print, date+' - all done!'

end
