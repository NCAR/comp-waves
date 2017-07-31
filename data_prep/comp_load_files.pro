;
;PURPOSE: Loads suitable CoMP files and interpolates. Called by wave_tracking.pro
;
;OPTIONAL INPUTS: /debug - for debugging
;                 lim  - sets number of missing data files allowed to interpolate over.
;                        Value is set in wave_tracking- default is 2
;                        code cannot deal with more than this yet (probably sensible anyway!) 
;
;OUTPUTS: cube_i - intensity cube
;         cube_v - velocity cube
;         cube_w - doppler width cube
;
;
;CALLS: file_list.pro,headfits.pro,sxpar.pro,readfits.pro
;
;HISTORY: Written by R J Morton 05/2015 using parts of wave_tracking.pro by S. Tomczyk
;         Made data interpolation routine generic - RJM 06/2015
;         Created index & temp_index to keep together useful quantities from wave tracking RJM 04/2016
;         Fixed bug for files with no missing frames - RJM/C. Bethge 06/2016

PRO comp_load_files,cube_i,cube_v,cube_w,index,debug=debug,lim=lim,cadence=cadence

COMMON wt_var,inpath,nwlst,rsun_mm, temp_index,ccthresh,maxiter,num_files_proc     ;start_file,norm_cadence,nx,ny,nt,num_files_proc,xscale

;flist   = file_list(inpath,'*.comp.1074.dynamics.'+nwlst+'.fts.gz')
flist   = file_list(inpath,'*.comp.1074.dynamics.fts.gz')
if temp_index.start_file ne 0 then begin
 flist = flist[temp_index.start_file:*]
endif
fnumber = n_elements(flist)
cadence = intarr(fnumber-1)
hdr     = headfits(flist[0],ext=1)
nx      = sxpar(hdr,'NAXIS1')
ny      = sxpar(hdr,'NAXIS2')

index=fitshead2struct(headfits(flist[0]) )


;================================================================
;compute cadence
;================================================================
FOR ii=0,fnumber-2 do begin
 hdr       = headfits(flist[ii])
 date_obs  = sxpar(hdr,'DATE-OBS')
 time_obs  = sxpar(hdr,'TIME-OBS')
 IF n_elements(time) EQ 0 THEN time=time_obs ELSE time=[temporary(time),time_obs]
 new_time  = anytim2tai(date_obs+' '+time_obs) 
 IF ii GT 0 THEN cadence[ii-1] = fix(new_time-old_time)
 old_time = new_time
  
ENDFOR


IF keyword_set(debug) THEN plot,cadence,xtitle='Frame No.',ytitle='Cadence (s)'
norm_cadence = median(cadence)

;================================================================
; look for the start of the normal cadence in case the sequence 
; starts with files with an abnormal cadence
; RJM updated section - MAY2015
;================================================================
seq_start = 0
seq_end   = fnumber-1

in=where(cadence EQ norm_cadence)
seq_start=in[0]


;================================================================
; interpolate over gaps of one missing file
; (happens quite often), otherwise don't bother
; RJM updated section - MAY2015
;================================================================

IF n_elements(lim) EQ 0 THEN lim=2

;find locations where cadence is greater than the normal cadence
;and locations where the jump is larger than the limit x normal cadence
in=where(cadence[seq_start:fnumber-(lim+1)] GT 1.1*norm_cadence)
in2=where(cadence[seq_start:fnumber-(lim+1)] GT (lim)*norm_cadence)

;Determine final image in series
IF in2[0] EQ -1 THEN seq_end=fnumber-(lim+1) ELSE seq_end=in2[0]

;Determine positions of missing frames
itrplt_pos=in[where(in lt seq_end)]
IF where(in[0] lt seq_end) EQ -1 THEN BEGIN 
   nmiss=0 
   itrplt_counter=0
ENDIF ELSE IF itrplt_pos[0] GE 0 THEN nmiss=n_elements(itrplt_pos) ELSE nmiss=0 ; number of missing frames

IF nmiss GT 0 THEN itrplt_counter=n_elements(itrplt_pos) ELSE itrplt_counter=0

;Works out which positions need filling
;including the extension of the array for previously 
;interpolated frames
IF nmiss GT 0 THEN BEGIN
      k=0
      itrplt_pos=0
      itrplt_val=0
      FOR i=0,itrplt_counter-1 DO BEGIN
          val=(cadence[seq_start:fnumber-(lim+1)])[in[i]]/norm_cadence
          
          CASE val OF
               2: BEGIN
                    itrplt_pos=[itrplt_pos,in[i]+k]
                    k=k+1
                    itrplt_val=[itrplt_val,val]
                  END   
               3: BEGIN
                    itrplt_pos=[itrplt_pos,in[i]+k,in[i]+1+k]
      	      k=k+2
                    itrplt_val=[itrplt_val,val,val]
                  END
               ELSE: message,"cannot yet interpolate over >4 frames" 
          ENDCASE
          
      ENDFOR
      itrplt_pos=itrplt_pos[1:itrplt_counter]
      itrplt_val=itrplt_val[1:itrplt_counter]
      print,'Missing frames',itrplt_pos
ENDIF 


;================================================================
; cut out contiguous file list and prepare
; the interpolated data cube
; RJM updated section - MAY2015
;================================================================
print,'Start sequence',seq_start,' End sequence',seq_end

;Cuts down file list to those values being loaded in
flist = flist[seq_start:seq_end]
IF nmiss GT 0 THEN itrplt_pos = itrplt_pos+1 ;+seq_start
nt = seq_end-seq_start+1+itrplt_counter

cube_i = fltarr(nx,ny,nt)
cube_v = fltarr(nx,ny,nt)
cube_w = fltarr(nx,ny,nt)


xscale = rsun_mm*index.CDELT1/index.RSUN ;Mm
time_test = strarr(nt)+'ZERO' ; (in)sanity check

;================================================================
; fill cubes with data and interpolate if necessary
; RJM updated section - MAY2015
; Should now be able to interpolate over any number of frame gaps…
;…IF DESIRED - 'Authenticity' decreases with increasing interpolation gaps
;================================================================

;First read in files and put in cube in their place, leaving gaps for missing frames
k=0
add=0
FOR ii=0,nt-1 DO BEGIN
   
    IF ii NE (itrplt_pos[k]) THEN BEGIN
        cube_i[*,*,ii] = readfits(flist[ii-add],ext=1,/silent)
        cube_v[*,*,ii] = readfits(flist[ii-add],ext=3,/silent)
        cube_w[*,*,ii] = readfits(flist[ii-add],ext=4,/silent)
        time_test[ii]  = sxpar(headfits(flist[ii-add]),'TIME-OBS')
    ENDIF ELSE add=add+1

    IF k LT itrplt_counter-1 THEN IF ii EQ itrplt_pos[k] THEN k=k+1

ENDFOR

;Fill in missing frames
ii=0
IF itrplt_counter GT 0 THEN BEGIN
    WHILE ii LT itrplt_counter DO BEGIN
         tt=itrplt_pos[ii]
         tval=itrplt_val[ii]
         z=findgen(tval)/tval
         z=z[1:tval-1]

         dumi=[[[cube_i[0:nx-1,0:ny-1,tt-1]]],[[cube_i[0:nx-1,0:ny-1,tt+tval-1]]]]
         dumv=[[[cube_v[0:nx-1,0:ny-1,tt-1]]],[[cube_v[0:nx-1,0:ny-1,tt+tval-1]]]]
         dumw=[[[cube_w[0:nx-1,0:ny-1,tt-1]]],[[cube_w[0:nx-1,0:ny-1,tt+tval-1]]]]
         cube_i[0:nx-1,0:ny-1,tt:tt+tval-2] = interpolate(dumi,findgen(nx),findgen(ny),z,/grid)
         cube_v[0:nx-1,0:ny-1,tt:tt+tval-2] = interpolate(dumv,findgen(nx),findgen(ny),z,/grid)
         cube_w[0:nx-1,0:ny-1,tt:tt+tval-2] = interpolate(dumw,findgen(nx),findgen(ny),z,/grid)
         ii+=tval-1
         
    ENDWHILE
  
ENDIF


; manually exclude bad data
IF num_files_proc NE 0 THEN BEGIN
  IF num_files_proc GT nt THEN num_files_proc=nt
  cube_i = reform(cube_i[*,*,0:num_files_proc-1])
  cube_v = reform(cube_v[*,*,0:num_files_proc-1])
  cube_w = reform(cube_w[*,*,0:num_files_proc-1])
  nt=num_files_proc
ENDIF


index=create_struct(temporary(index), 'NAXIS1', nx, 'NAXIS2', ny, 'NFRAMES', nt, $ 
                    'TIME_D$OBS_LIST',time_test[0:nt-1], 'XSCALE', xscale, 'NORM_CADENCE',  norm_cadence, temp_index)



END