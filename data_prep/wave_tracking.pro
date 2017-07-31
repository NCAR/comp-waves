;  PURPOSE: Procedure to determine the wave propagation angle from the cross
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
;              lim      - limit to number of missing frames in a row, default is 2 (cannot do more than this)
;              dosob    - uses sobel filtered image for alignment - better & faster results
;              splt_cube - save cubes as separate files
;              npt_init - choose path length for phase speed measurement (default 25) - choose odd value
;
; CALLS: find_man_date.pro, do_apod.pro, comp_load_files.pro, wave_angle_calc.pro, compute_speed.pro, dejitter.pro,
;        plot_wave_maps.pro
;
;
;
;
;
;HISTORY - written by S. Tomczyk & C. Bethge
;          Modularised by R J Morton 2015
;          Made data interpolation routine generic - RJM 06/2015
;          Reduced computation time for wave_angle_calculation - RJM 06/2015
;          Added apodisation window for wave_angle_calculation, periodic nature of FFT means non-zero endpoints
;          leads to spectral leakage, can affect coherence estimates. Apodisation minimises leakage.
;          Not clear which apodistaion window is best to use for CC just yet - RJM 06/2015
;          Added option to use Sobel filtered data for alignment - RJM 11/2015
;          Created temp_index & index to keep together useful quantities from wave tracking - RJM 04/2016
;          Redefinition of mask values using coherence measure form wave angles - RJM + CB 06/2016
;


pro wave_tracking_v2, date, nwl=nwl,                           $
                         own_data=own_data,                 $
                         compute_speeds=compute_speeds,     $
                         cross_corr=cross_corr,             $
                         plot_maps=plot_maps,               $
                         save_coh = save_coh,               $
                         name_addon = name_addon,           $
                         wr_speeds = wr_speeds,             $
                         choose_corr_box = choose_corr_box, $
                         debug = debug, altangle=altangle,  $
                         freq0=freq0, $
                         lim=lim, dosob=dosob,splt_cube=splt_cube,npt_init=npt_init


COMMON wt_var,inpath,nwlst,rsun_mm, temp_index,ccthresh,maxiter,num_files_proc     ;start_file,norm_cadence,nx,ny,nt,xscale
COMMON wt_var_ang,mask,dx,nbox,nsmooth,limit,npix,fwidth,freq_filt,good

;================================================================
;*** Here the details of the wave tracking need to be defined ***
;================================================================

temp_index=create_struct('NWL', 0, 'START_FILE', 0, $
                       'LOWER_R', 0, 'UPPER_R', 0, $
                       'MAXOFFSET', 0.0, 'CC_COORD', fltarr(2,2), 'RMS_CC', fltarr(2))

; use 3 pt data by default
if n_elements(nwl) eq 0 then temp_index.nwl = 3 $
ELSE temp_index.nwl=nwl
nwlst = strcompress(string(temp_index.nwl),/rem)
ans   = ' '

rsun_mm = 695.7 ;solar radius in Mm

if keyword_set(wr_speeds) then wr_speeds = 1 else wr_speeds = 0

; paths where the CoMP data is located and where the .sav-files should be written to
inpath  = 'analysis/CoMP/'+strmid(date,8,4,/reverse)+'/'+strmid(date,3,2,/reverse)+'/'+strmid(date,1,2,/reverse)+'/'
outpath = 'analysis/CoMP/wave_tracking_output/'+date+'/'


if file_test(outpath,/DIR) eq 0 then spawn, 'mkdir '+outpath
if keyword_set(plot_maps) then begin
 if file_test(outpath+'plots/',/DIR) eq 0 then spawn, 'mkdir '+outpath+'plots/'
endif

;addon to the names of the .sav-files
if n_elements(name_addon) eq 0 then name_addon = ''else name_addon = name_addon+'_'


; values for finding the wave propagation angles
dx      = 20      ;
nbox    = 2*dx+1  ;box size in pixels
nsmooth = 15      ;smoothing width for cross spectra (nominally 11)
limit   = 0.5     ;set coherence limit for acceptance
npix    = 10      ;minimum number of coherent pixels
fwidth = 0.0015

IF n_elements(freq0) LT 1 THEN freq0  = 0.0035
freq_filt=freq0


case freq0 of
 0.0015: freqst = '1.5mHz'
 0.0035: freqst = '3.5mHz'
 0.0055: freqst = '5.5mHz'
 else: begin
         freqst = 'manual_freq_setting'
       end
endcase


;===============================================
;Processing values for default and certain dates
;===============================================

find_man_date,date ;,date,x1,x2,y1,y2,start_file,num_files_proc,lower_r,upper_r,ccthresh,maxiter
box = temp_index.cc_coord




;====== if keyword own_data is set  =================================
; skip the preparation of the data if someone provides his own
; velocity data cube (cube_v) - in this case, a few additional 
; variables need to be defined here: 
; norm_cadence = the cadence of cube_v in seconds
; arcsec_per_pixel = usually 4.46 for CoMP data
if n_elements(own_data) ne 0 then begin
 print,'Restoring user data'
 restore, own_data
 s=size(cube_v)
 nx=s[1] 
 ny=s[2]
 nt=s[3]
 
 norm_cadence     = 30.
 arcsec_per_pixel = 4.46
 xscale              = arcsec_per_pixel*0.725   
                    ;pixels in Mm - assume 725 km/arcsec (1 AU distance). 
                    ;correct this value if necessary -
                    ;important for correct phase speeds!!
 goto, external_data
endif


;================================================================
;***      End of defining the wave tracking parameters        ***
;================================================================


;================================================================
;READ IN FILES
;================================================================
comp_load_files,cube_i,cube_v,cube_w,index,debug=debug,lim=lim,cadence=cadence



external_data:

;================================================================
; cross-correlate the data and align it
;================================================================
if keyword_set(cross_corr) then begin
  if keyword_set(debug) then ccdebug = 1

  if keyword_set(choose_corr_box) then begin
     IF keyword_set(dosob) THEN im=sobel(cube_i[*,*,0])<20 ELSE im=cube_i[*,*,0] >0 <25
     ONCE_MORE:
     x2 = -1
     y2 = -1
     window, /free, xs=620, ys=620
     win_nr = !d.window
     aia_lct, wave=193, /load
     tvscl, reform(im)
     message, /cont, 'Choose bottom left corner'
     cursor, x1, y1, /down, /device    
     message, /cont, 'Choose top right corner'

     IF !MOUSE.BUTTON EQ 1 THEN BEGIN 
        CURSOR, do_nothing_x, do_nothing_y, /UP, /WAIT
        !MOUSE.BUTTON = 0  
     ENDIF
     while !Mouse.Button NE 1 do begin
        cursor, tmp_x2, tmp_y2, /change, /device
        tvscl, reform(im)
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
     index.CC_COORD=[[x1,y1],[x2,y2]]
  endif

  print,'Beginning Alignment'
  dejitter,cube_i,cube_v,cube_w,index,ccdebug=ccdebug,dosob=dosob
  
endif


print,'Stop to check alignment'
print,'Type .c to continue'
stop


message, /cont, 'Total # of files: '+strcompress(string(index.NFRAMES,format='(F5.0)'),/rem)

;================================================================
;FFT of velocity data cube
;RJM - 29 May 2015 - rearranged masking and FFT
; used mask to reduce number of time-series FFT'd by
; removing the empty time-series ~75% less operations required!
;================================================================

print,'Begin FFT of data cube

; Create mask - mask out pixels not to invert
nx=index.NAXIS1
ny=index.NAXIS2
mask=intarr(nx,ny)
x=rebin(findgen(ny)-(nx/2.)+0.5,nx,ny)
y=transpose(x)
r=sqrt(x^2+y^2)
good=where(r ge index.LOWER_R+index.MAXOFFSET and r le index.UPPER_R)
mask(good)=1



; open windows
if keyword_set(debug) then begin
 loadct, 39, /silent
 f=8
 window,1,xs=nbox*f,ys=nbox*f,xpos=0,ypos=0
 window,4,xs=nx,ys=ny,retain=2,xpos=400+nx+20,ypos=0, title=date
 window,5,xs=ny,ys=ny,retain=2,xpos=400,ypos=ny+50
 window,6,xs=512,ys=512,xpos=400,ypos=800
 device,decomposed=0
endif





;================================================================
;START THE CALCULATION OF WAVE ANGLE
;================================================================

nspec = index.NFRAMES/2
spec=complexarr(nx,ny,nspec)
h=do_apod(index.NFRAMES,cpg) ; apodisation window RJM 06/2015

FOR iy=0,ny-1 DO $
FOR ix=0,nx-1 DO BEGIN

  IF mask[ix,iy] EQ 1 THEN BEGIN ;RJM 29052015
	  d=reform(cube_v[ix,iy,*]-mean(cube_v[ix,iy,*])) 
  	sp=fft(d*h)                ;multiply by apodisation window RJM 06/2015
    spec[ix,iy,*]=sp[0:nspec-1]
  ENDIF
ENDFOR

wave_angle_calc,spec,index,wave_angle,angle_error,coh_measure,save_coh=save_coh,debug=debug


;update mask based on coherence results from wave angle calculation
coh_qual = (coh_measure[*,*,0]-coh_measure[*,*,1])/(coh_measure[*,*,0]+coh_measure[*,*,1])
coh_qual[where(finite(coh_qual) EQ 0)] = 0.
mask[where(coh_qual LE 0)]=0.


index=create_struct(temporary(index), 'MASK', mask,'COH_MEAS',coh_qual)

IF keyword_set(save_coh) THEN save, coh_measure, file=outpath+'coherence_measure_'+name_addon+date+'_'+nwlst+$
                                    '_'+freqst+'.sav'

save,wave_angle,angle_error,index,file=outpath+'wave_angle_'+name_addon+date+'_'+nwlst+'_'+freqst+'.sav'

IF keyword_set(splt_cube) THEN BEGIN
    save, cube_i, index, filename=outpath+'cube_i_'+date+name_addon+'.sav'
    save, cube_v, index, filename=outpath+'cube_v_'+date+name_addon+'.sav'
    save, cube_w, index, filename=outpath+'cube_w_'+date+name_addon+'.sav'
ENDIF ELSE $
save, cube_i, cube_v, cube_w, index, filename=outpath+'cube_ivw_'+date+name_addon+'.sav'


IF keyword_set(compute_speeds) THEN BEGIN
 print,'Calculating propagation velocities'
 space_time_run, date, cube_v, wave_angle, index, outpath, nwlst, freqst, name_addon, $
                        debug=debug, wr_speeds=wr_speeds,npt_init=npt_init
ENDIF


IF keyword_set(plot_maps) THEN BEGIN
   IF  keyword_set(save_coh) THEN BEGIN
       plot_wave_maps, date, nwl, freqst, outpath, name_addon, /cohm
   ENDIF ELSE BEGIN 
       plot_wave_maps, date, nwl, freqst, outpath, name_addon
   ENDELSE 
ENDIF

print, date+' - all done!'

END
