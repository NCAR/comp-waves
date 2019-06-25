;PURPOSE: Calculates the angle of wave propagation from a data cube.
;         Current version is designed for CoMP velocity data.
;
;
;
;INPUTS: spec - cube of FFT'd data that runs from f=nt/2 to f=0 (i.e., half cube)
;
;OPTIONAL INPUTS: /save_coh - saves coherence values
;        
;OUTPUTS: wave_angle - array of wave angles
;         angle_error - error on wave angles
;         coh_measure - measure of coherence, only created if /save_coh is set
;
;
;NOTES: Array defined call mask. If STEREO idl files are installed then routine can get confused
;       and tries calling the STEREO routine. Maybe should change the name of mask array RJM 06/2015
;
;HISTORY: Created by ???
;         Modularised by R J Morton 2015
;         Created index to keep together useful quantities from wave tracking RJM 04/2016         
;
;TO DO: Is the current method of calculating coherence correct? Should summing be used instead of smoothing? - RJM
;


PRO wave_angle_calc_v02,spec,index,wave_angle,angle_error,coh_measure,save_coh=save_coh,debug=debug

COMMON wt_var,inpath,nwlst,rsun_mm, temp_index,ccthresh,maxiter,num_files_proc 
COMMON wt_var_ang,mask,dx,nbox,nsmooth,limit,npix,fwidth,freq_filt,good

nx=index.NAXIS1
ny=index.NAXIS2
nt=index.NFRAMES

wave_angle=fltarr(nx,ny)
angle_error=fltarr(nx,ny)
angle_error2=fltarr(nx,ny) ; For testing RJM
coh_measure = fltarr(nx,ny,2)
f=8

xbox=rebin(findgen(nbox)-dx,nbox,nbox)  ;coordinates of points in box
ybox=transpose(xbox)

;Set up filter
nspec=nt/2
freq   = findgen(nspec)/(float(nspec*2)*index.NORM_CADENCE)
filter = exp( -(freq-freq_filt)^2/fwidth^2 )
filter(0) = 0.                      ; set dc to zero
filter(where (freq lt .001)) = 0.   ; set low frequencies to zero
filter=filter/total(filter)



if not keyword_set(debug) then begin
 pixcounter = 0L
 old_perc_proc = 0
 tot_pix = long(n_elements(good))
endif

conjspec=conj(spec)
lar_g=real_part(smooth(spec*conjspec,[1,1,nsmooth]))

message, /cont,'Starting processing...'
FOR iy=0,ny-1 DO FOR ix=0,nx-1 DO IF mask(ix,iy) EQ 1 THEN BEGIN

  IF NOT keyword_set(debug) THEN BEGIN
     perc_proc = round((float(pixcounter)/float(tot_pix))*100.)
     pixcounter = pixcounter+1
     IF perc_proc GT old_perc_proc THEN BEGIN
            IF strmid(getenv('TERM'),0,5) EQ 'xterm' THEN BEGIN
     		progress_comp, perc_proc, 100, msg='Computing wave propagation angles'
    	    ENDIF ELSE message, /cont, strcompress(string(perc_proc,format='(I3)'),/rem)+'% processed.'
    	    old_perc_proc = perc_proc
     ENDIF
  ENDIF
 
  ;================================================================
  ;  compute coherence for contiguous pixels with coherence greater than limit
  ;================================================================

  coh=fltarr(nbox,nbox)           ;initialize coherence array to zero
  coh_mask=fltarr(nbox,nbox)      ;initialize mask of good points in coherence
  coh(nbox/2,nbox/2)=1.
  coh_mask(nbox/2,nbox/2)=1.

  spec1=reform(spec[ix,iy,*])
  ;spec1c=reform(conjspec[ix,iy,*])
  ;g1=real_part(spec1*spec1c)
  ;g1=smooth(g1,nsmooth)
  g1=reform(lar_g[ix,iy,*])

  count=1
  ncoh=1
  WHILE count GT 0 AND ncoh LT 100 DO BEGIN
        cont=where(coh_mask eq 1.)
        cnew=[cont-1,cont+1,cont+nbox,cont-nbox]  ;test pixels adjacent to contiguous pixels

        new_mask=fltarr(nbox,nbox)   ;identify new contiguous pixels
        new_mask(cnew)=1.
        new_mask=new_mask*(1.-coh_mask)

        new=where(new_mask gt 0.,new_count)
       
        count=0
        FOR ic=0,new_count-1 DO BEGIN
            icx=new(ic) mod nbox
            icy=new(ic)/nbox

            ii=ix-nbox/2+icx
            jj=iy-nbox/2+icy

            IF ii LT 0 OR ii GT nx-1 OR jj LT 0 OR jj GT ny-1 THEN coh(icx,icy)=0. ELSE BEGIN
        	  IF coh(icx,icy) EQ 0. THEN BEGIN
          	         ;spec2=reform(spec(ii,jj,*))
                     spec2c=reform(conjspec(ii,jj,*))
                     cspec=spec1*spec2c
                     cspec=smooth(cspec,nsmooth)
                     g2=reform(lar_g[ii,jj,*])
                     ;g2=real_part(spec2*spec2c)
                     ;g2=smooth(g2,nsmooth)
                                          
                     coh(icx,icy)=mask(ii,jj)*total( filter*( abs(cspec)/sqrt(g1*g2) ) ) 
                     ;see McIntosh et al. 2008. Sol. Phys.
                  ENDIF

                  IF coh(icx,icy) GT limit THEN BEGIN
          	     coh_mask(icx,icy)=1.
          	     count=count+1
                  ENDIF
            ENDELSE
         ENDFOR

         good=where(coh_mask eq 1.,ncoh)

   ENDWHILE

   IF keyword_set(debug) THEN BEGIN
    wset,1
    tv,bytscl(rebin(coh*coh_mask,f*nbox,f*nbox,/sample),0.,1.)
   ENDIF   

   ;IF (keyword_set(save_coh) and (ncoh gt npix)) THEN BEGIN
   IF (ncoh gt npix) THEN BEGIN
     island = coh*coh_mask
     ellipsepts = fit_ellipse(where(island ne 0.),axes=axes,xsize=nbox,ysize=nbox)
     coh_measure[ix,iy,*] = axes
     if ((size(ellipsepts))[0] gt 0) and keyword_set(debug) then plots, f*ellipsepts, color=254, /device
   ENDIF

   ;================================================================ 
   ;  find angle of coherence island if enough good points
   ;================================================================

   IF ncoh GT npix THEN BEGIN

 
       weight=coh(good)^2
       theta=min_dist_fit(xbox(good),ybox(good),error=error) ; analytical solution to the minimum distance problem
       ;theta=mle_angle(xbox(good),ybox(good)) 
       ;error=0.

       IF keyword_set(altangle) THEN BEGIN

           xe=xbox(good) & ye=ybox(good)
           sxy=total((xe-mean(xe))*(ye-mean(ye)))

           IF round(sxy) NE 0. THEN BEGIN
    
                 sixlin,xbox(good),ybox(good),a,siga,b,sigb
                 theta=[atan(b[2])/(!dtor)]
                 error=[atan(sigb[2])/(!dtor)]
           ENDIF
       ENDIF
        
       wave_angle(ix,iy)=theta
       angle_error(ix,iy)=error


       
       IF keyword_set(debug) THEN BEGIN
          xfit=findgen(2*dx+1)-dx
          coef=[0.,tan(theta*!dtor)]
          yfit=poly(xfit,coef)  ;C0 + C1*x + C2*x^2

          wset,6
          plot,xbox(good),ybox(good),psym=6,xr=[-dx,dx],yr=[-dx,dx]
          oplot,xfit,yfit
          wset,4
          tv,bytscl(wave_angle,-90,90)
          wset,5
          tv,bytscl(angle_error,0,10)
          wait,0.005
          print,theta,error
       ENDIF

  ENDIF
  
ENDIF


END
