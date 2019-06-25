;PURPOSE - Create a time-distance diagram along an arc in CoMP data. 
;
;INPUTS - indir - directory of data cubes. 
;                 Requires data to be saved as intensity.fcube,velocity.fcube,width.fcube
;                 Also requires ROI to be saved as intensity_loop.fcube,…
;                 Also requires file called smallwindowinfo.sav which contains: x1box - x coordinates of cutout    
;									     y1box - y coordinates of cutout
;  									     initime - initial time of observation
;									       
;OPTIONAL INPUTS: cutsize - half length of cut perpendicular to arcade,default=10 
;                 jump - can skip samples along arc if required,default=1
;                 upix - pixels in unsharp mask, default=5
;
;
;
;HISTORY: Created by J. Threfall (2012)
;         Tidied & edited basic functionality by R Morton (2014)



FUNCTION pick_path,p

;At present works only for loops on east limb
;Have to very accurate with selection else you get a wonky arc
  FOR pc=0,5 do BEGIN
    pick,x,y
    p[*,pc]=[x,y]
    plots,[x,x],[y,y],psym=1,thick=3
  ENDFOR
  out=[[p[0,0]+10,p[1,0]-10], [p[0,*],p[1,*]], [p[0,5]+10,p[1,5]+10]]
  plots,[out[0,0],out[0,0]],[out[1,0],out[1,0]],psym=1,thick=3
  plots,[out[0,7],out[0,7]],[out[1,7],out[1,7]],psym=1,thick=3

return,out

END

FUNCTION define_path,pextra=pextra,x0c=x0c,y0c=y0c,dx=dx,dy=dy
     ;First off, pick some global coordinates along the thing you wish to follow.
     ; you can use cursor, x, y, /data, or copy in some coordinates from another instrument (like I did here)
     ; p IS THE ARC COORDINATES:
     ; pextra goes beyond this range in order to get the shape of the arc right after ArcSample
     ; pextra poinst should be placed around +/-[10,10] pixels away depending on path

     ;coordinates for 20120410 p=[[47.25,44.7],[33,53.4],[21.5,67],[19.5,87],[28,105.5],[44.8,116.25]]
     p=[[82.25,61.6],[63.82,77.6],[56.83,93.93],[58.69,109.33],[67.55,121.46],[78.19,129.5]]
     ;pextra=[[68,33], [p[0,*],p[1,*]], [65,125]]
     pextra=[[105,57], [p[0,*],p[1,*]], [102,135]]
     IF KEYWORD_SET(pickpath) THEN pextra=pick_path( p )

     ;-SPLINE- 
     ; we now have our coordinates through which we want a nice track - I will fit a spline so that we 
     ; have a track thats continuous, can bend, and follows the thing we want it too. 
     psize=size( p ) & pextrasize=size(pextra)
     xvals=pextra[0,psize[2]-1+2]+findgen(pextra[0,0]-pextra[0,psize[2]-1+2]) 
     spline_p, pextra[0,*], pextra[1,*], yvals, xvals

     ;-RESAMPLING-
     ; the problem with this is that the coordinates are equally spaced in y (or x), but not along the arc
     ; solution - arcsample program in the coyote IDL library resamples an arc to contain a specific number of points.
     ; here now follows a complicated bit - I will repeatedly run the arcsampling code, and increase the number of points, 
     ; until the distance between each point is LTE the CoMP pixel width.
     aseed=100 ; initial guess for the number of points
     result=0
    
     REPEAT BEGIN
       ;Arcsample creates closed loop
       ;Need to select points between arc endpoints 
       ArcSample, xvals, yvals, yinterp, xinterp, points=aseed
       xpoints=where(xinterp le pextra[0,1] and yinterp lt pextra[1,pextrasize[2]-2] and yinterp ge pextra[1,1])	    	    ; TAR arc case
  
       x0c=xinterp[xpoints] & y0c=yinterp[xpoints]
       dx=shift(x0c,-1)-x0c & dy=shift(y0c,-1)-y0c
       rad=sqrt(dx*dx+dy*dy) ; point separation distance
 
       IF rad[0] lt 1 THEN result=1     ; when less than 1 pixel width (=1!), exit
       aseed=aseed+1    	    	    ; otherwise increase number of sample points
     ENDREP UNTIL result
     dx=dx[0] & dy=dy[0]

return,p
END




PRO jt_compspline,indir=indir,pickpath=pickpath,cutsize=cutsize,jump=jump,upix=upix

;##############################
;FILE LOCATIONS
;##############################

;location of CoMP data folder
IF NOT keyword_set(indir) THEN Cdataloc='data/CoMP/wave_tracking_output/20120411/' ELSE Cdataloc=indir   

;location of saved spline data 
rottdloc=Cdataloc+'SPLINEcuts/'

;make directory if it doesn't exist	
res=findfile(rottdloc,count=count)
IF count EQ 0 THEN file_mkdir,rottdloc	

; names of the cubes for both full disk (FD) and small (S) subwindow:
CoMPFDcube=strarr(3)     	    	    	    
CoMPFDcube=Cdataloc+['Intensity.fcube','Width.fcube','Velocity.fcube']
CoMPFDcubeI=Cdataloc+'times.sav'

CoMPScube=strarr(4)	
CoMPScube[0:2]=Cdataloc+['intensity_loop.fcube','width_loop.fcube','velocity_loop.fcube']	   
;CoMPScube[3]=Cdataloc+'chris_SW_EI_DJ.fcube'	    ; velocity

;save file containing the x1, x2, y1, y2 coords of small window, the observation times
CinfoS=Cdataloc+'smallwindowinfo.sav'
restore, CinfoS, /verbose


;############################################################
;Main routine
;############################################################
;Set plot windows
WINDOW, 0, ysize=500, xsize=400
wset,0
WINDOW, 1, ysize=500, xsize=400 
wset,1

  
IF n_elements(cutsize) EQ 0 THEN cutsize=10 ; Defines half length of cut perpendicular to arc
IF n_elements(jump) EQ 0 THEN jump=1     ; jump tells you the skip between each sample (increasing above one increases speed, at a cost of info).
IF n_elements(upix) EQ 0 THEN Upix=5  	    	    	    	

LOOPTRACKER=1	    	; turns on the actual spline fitting bit (if =1)

IF LOOPTRACKER THEN BEGIN
  !p.background=255

 FOR q=0,2 DO BEGIN
  ;setting up different directories for each variable, in case one overwrites another. 
  CASE q OF
   0: BEGIN
       compdat='I'
       subdir=rottdloc+'intensity/'
       res=findfile(subdir,count=count)
       IF count EQ 0 THEN file_mkdir,subdir	
      END
   1: BEGIN
       compdat='LW'
       subdir=rottdloc+'linewidth/'
       res=findfile(subdir,count=count)
       IF count EQ 0 THEN file_mkdir,subdir
      END
   2: BEGIN
       compdat='V'
       subdir=rottdloc+'velocity/'
       res=findfile(subdir,count=count)
       IF count EQ 0 THEN file_mkdir,subdir
      END
   3: BEGIN
       compdat='EI'
       subdir=rottdloc+'enhancedintensity/'
       res=findfile(subdir,count=count)
       IF count EQ 0 THEN file_mkdir,subdir
      END
   ELSE: print, 'CoMP data selection failure'
  ENDCASE

  rotslice=subdir+'arc1.sav'	    ; name of the actual save file containing the spline data
  
  
  asize=size(lp_read(CoMPScube[q])) ; size of the CoMP window we're using
  
  
  wset,0
  loadct, 0,/silent
  pih, alog10(lp_get(CoMPScube[q],0)),$
  title='CoMP cutout, '+compdat+', '+initime, col=0;, $
  
  ;Only need to run once
  IF q EQ 0 THEN p=define_path(pextra=pextra,x0c=x0c,y0c=y0c,dx=dx,dy=dy)

  psize=size(p)
  pextrasize=size(pextra)
  FOR i=0,pextrasize[2]-1 DO BEGIN
   plots, [pextra[0,i],pextra[0,i]], [pextra[1,i],pextra[1,i]], thick=3, psym=7, symsize=2
   xyouts, pextra[0,i]+2,pextra[1,i]+2,strcompress(i-1), charthick=2, charsize=2
  ENDFOR 
  FANcoords=pextra
  plots, x0c, y0c, thick=2	    	    ; resulting coords, spaced with the width of one pixel.
  stop

  ;############################################################
  ;-TAKING PERP CUTS AT EACH POSITION
  ;############################################################

  dat=lp_read(CoMPScube[q])	   
  ds=size(dat)
  cco=[ds[1]/2,ds[2]/2]	    	    	; coords of centre
  cutx=[cco[0]-cutsize,cco[0]+cutsize]	; coordinates to cut from and too
  nsamples=(n_elements(y0c)-jump)/jump
  slices=dblarr(2*cutsize-1,6,ds[3]+1,nsamples)   ;define the place we will store the data
  tot=0     	    	    	    	; count, progress bar
     
  ; LOOP OVER EACH POSITION ALONG THE ARC AND SAMPLE  
  FOR j=0,n_elements(y0c)-jump-1,jump DO BEGIN 
   ;work out the angle the arc makes with the vertical
   IF x0c[j+jump] lt x0c[j] THEN BEGIN    	    	    
      ; if in normal quad, then angle=atan(dx/dy)
      IF y0c[j+jump] gt y0c[j] THEN theta=atan(ABS(DOUBLE(x0c[j+jump]-x0c[j]))/ABS(DOUBLE(y0c[j+jump]-y0c[j]))) $
      ELSE theta=!dpi-atan(ABS(DOUBLE(x0c[j+jump]-x0c[j]))/ABS(DOUBLE(y0c[j+jump]-y0c[j]))); if in other quad, angle=pi-atan(dx/dy)

   ENDIF ELSE BEGIN
     ; third quad, angle=-atan(dx/dy)
     IF y0c[j+jump] gt y0c[j] THEN theta=-atan(ABS(DOUBLE(x0c[j+jump]-x0c[j]))/ABS(DOUBLE(y0c[j+jump]-y0c[j]))) $   
     ELSE theta=!dpi+atan(ABS(DOUBLE(x0c[j+jump]-x0c[j]))/ABS(DOUBLE(y0c[j+jump]-y0c[j]))); fourth quad, angle=pi+atan(dx/dy)
   ENDELSE

   ; Rotate the entire cube so that the cut doesn't have to move.
   ; guarantees to be cutting exactly perpendicular to the arc
   datr = CUBEROT( dat, theta*!radeg, x0c[j],y0c[j],missing=0)
   
   wset,1
   loadct, 0, /silent
   pih, datr[*,*,0], title=compdat+string(tot,format='(", cut ",I3.3)'), col=0     ; shows the orientation and position of new cube
  
   loadct, 1,/silent
   ; arrow highlights where the spline fit is and current angle
   arrow, cco[0],cco[1]-30,cco[0],cco[1]-2, color=100, /data , thick=3	    	    ; arrows highlight location we are cutting
   arrow, cutx[0]-10,cco[1],cutx[0],cco[1], color=230, /data , thick=2
   arrow, cutx[1]+10,cco[1],cutx[1],cco[1], color=230, /data , thick=2
   xyouts, cutx[0]-30, cco[1], 'cut:', col=230 
   
  
   ;now begin sampling values - take the t=0 ones first (running diff needs extra step to get going)
   ; cutter subroutine simply samples the data grid, from x0 to x1, y0 to y1
   ; note y0 and y1 are the same:
   tempo=cutter(datr[*,*,0], cutx[0],cutx[1], cco[1], cco[1])
   ;temps=cutter(alog10(datr[*,*,0])-alog10(smooth(datr[*,*,0],Upix,/edge_truncate)), cutx[0],cutx[1], cco[1], cco[1])
   ;tempss=cutter(alog10(datr[*,*,0])-smooth(alog10(datr[*,*,0]),Upix,/edge_truncate, /NAN), cutx[0],cutx[1], cco[1], cco[1]) 
  
   
   nu=cutx[1]-cutx[0]-2
   
   slices[*,0,0,tot]=tempo[0:nu,1]	    ; original intensity variation
   ;slices[*,1,0,tot]=temps[0:nu,1]	    ; USM (alog(im-smooth(im))
   ;slices[*,2,0,tot]=0	    	            ; rdiff
   ;slices[*,3,0,tot]=tempss[0:nu,1]         ; newer USM (alog(im)-smooth(alog(im))

   ;now loop over time and perform the rest of the samples
   FOR jj=1,asize[3]-1 DO BEGIN
    tempo=cutter(datr[*,*,jj], cutx[0],cutx[1], cco[1], cco[1])
    ;temps=cutter(alog10(datr[*,*,jj])-alog10(smooth(datr[*,*,jj],Upix,/edge_truncate)), cutx[0],cutx[1], cco[1], cco[1])
    ;tempss=cutter(alog10(datr[*,*,jj])-smooth(alog10(datr[*,*,jj]),Upix,/edge_truncate,/NAN), cutx[0],cutx[1], cco[1], cco[1]) 
    ;tempr=cutter(datr[*,*,jj]-datr(*,*,jj-1), cutx[0],cutx[1], cco[1], cco[1])
      
    slices[*,0,jj,tot]=tempo[0:nu,1]	    	; original intensity variation
    ;slices[*,1,jj,tot]=temps[0:nu,1]	    	; subnet mask
    ;slices[*,2,jj,tot]=tempr[0:nu,1]	    	; rdiff
    ;slices[*,3,jj,tot]=tempss[0:nu,1]  	       	; new USM

   ENDFOR
  
   tot=tot+1
   ;print, 'tot=', tot, '/', nsamples
   wset,0
   loadct, 39, /silent
   plots, x0c[j],y0c[j], col=100, psym=7,thick=2     ; highlight current location on original track
  ENDFOR
  
  
  save, slices, x0c,y0c, times, FANcoords, filename=rotslice, /verbose, /compress
  Icube=subdir+'rotTDR'+'_Intenscube_'+compdat+string(jump,format='(I1)')+'.fcube'
  ;USMcube=subdir+'rotTDR'+'_USMcube_'+compdat+string(jump,format='(I1)')+'.fcube'
  ;USMcube2=subdir+'rotTDR'+'_USMcube2_'+compdat+string(jump,format='(I1)')+'.fcube'
  ;RDcube=subdir+'rotTDR'+'_RDcube_'+compdat+string(jump,format='(I1)')+'.fcube'


  lp_write, float(reform(slices[*,0,*,*])), Icube
  ;lp_write, float(reform(slices[*,1,*,*])), USMcube
  ;lp_write, float(reform(slices[*,2,*,*])), RDcube
  ;lp_write, float(reform(slices[*,3,*,*])), USMcube2

 ENDFOR
 
 ENDIF
 !p.background=0
END
