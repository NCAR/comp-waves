
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

FUNCTION compute_speed_new, vel, delx, delt, debug=debug, ret=ret,no_tr=no_tr


COMMON located_dat,located, nx, nt
COMMON threads_dat, threads

s=size(vel)
nt=s(1)
npt=s(2)
IF npt mod 2 EQ 0 THEN npt=npt-1 ;deals with even length tracks
icent=npt/2


atrous,vel,decomp=a,n_scales=1
vel=a[*,*,0] ;+smooth(a[*,*,1],3,/edge_truncate)



IF n_elements(algn_err) GT 0 THEN print,'Error on alignment is:', algn_err


numg=0.001
;READ,numg,PROMPT='Enter gradient: '

numt=5
;READ,numt,PROMPT='Enter minimum thread length: '


    IF not keyword_set(gauss) THEN BEGIN
      ;output is in common located_dat
      locate_things_min,data=vel[*,*],grad=numg
      
      ;print,'no Gaussian fit used'
    ENDIF ELSE BEGIN
      IF n_elements(errors) gt 0. THEN errorsi=errors[*,*]
   
      locate_things,data=vel[*,*],grad=numg,meas_size=meas_size, $
                 errors=errorsi,check=check,cut_chisq=cut_chisq
    
      IF n_elements(algn_err) GT 0 THEN located.errs=located.errs+algn_err   ; adds alignment error to error from Gaussian fits

      ;print,'Gaussian fit used'
    ENDELSE
    
    ;output is in common thread_dat
    follow_thread,min_tlen=numt,area=3
    threads_max=threads

    ;Do for minimum strips
    IF not keyword_set(gauss) THEN BEGIN
      ;output is in common located_dat
      locate_things_min,data=vel[*,*]*(-1.),grad=numg
      
      ;print,'no Gaussian fit used'
    ENDIF ELSE BEGIN
      IF n_elements(errors) gt 0. THEN errorsi=errors[*,*]
   
      locate_things,data=vel[*,*]*(-1.),grad=numg,meas_size=meas_size, $
                 errors=errorsi,check=check,cut_chisq=cut_chisq
    
      IF n_elements(algn_err) GT 0 THEN located.errs=located.errs+algn_err   ; adds alignment error to error from Gaussian fits

      ;print,'Gaussian fit used'
    ENDELSE
    
    ;output is in common thread_dat
    follow_thread,min_tlen=numt,area=3
    threads_mn=threads
   
    
    IF n_elements(threads_mn) gt 1 THEN threads=[threads_max,threads_mn[1:n_elements(threads_mn)-1]] $
    ELSE threads=threads_max


    sz=size(threads)
    len=0
    IF sz[1] GT 1 THEN BEGIN
     	 ph=fltarr(sz[1]-1,5)
     	 FOR i=1,sz[1]-1 DO BEGIN
        	dummy=threads[i].pos
        	terr=threads[i].err_pos
        	;if value missing set to same as last pixel and set
        	zer=where(dummy EQ 0.)
        	IF zer[0] NE -1 THEN BEGIN
           		FOR ii=0,n_elements(zer)-1 DO BEGIN
                    		dummy[zer[ii]]=dummy[zer[ii]-1]
                    		terr[zer[ii]]=1.
           		ENDFOR
        	ENDIF
        	in=where(dummy gt -1.)
        	;len=[len,n_elements(in)]
        	res=poly_fit(findgen(n_elements(in)),dummy[in],1,yfit=fit,measure_err=terr[in],sigma=sig,chisq=chi)
        	ph[i-1,0]=res[1] ;gradient
          ph[i-1,1]=sig[1] ;error
          ph[i-1,2]=chi/(n_elements(in)-2) ;reduced chi^2
          ph[i-1,3]=n_elements(in)         ;length of track
          ph[i-1,4]=in[0]
	     ENDFOR
    	 ;len=len[where(ph gt 0)]
       valf=1. ;Cutoff phase speed in Mm/s
       cutoff=1./(valf*delt/delx)

    	 IF keyword_set(ret) THEN ph=ph[where(ph[*,0] lt -1.*cutoff),*]  $
    	 ELSE ph=ph[where(ph[*,0] gt cutoff),*] 


       ;##############################################################
       ;bootstrapping of means
       ;##############################################################
       numbs=500
       vals=n_elements(ph[*,0])
      	
       IF vals GT 3 THEN BEGIN ;Arbitrary cut-off of minimum 3 values
      	      g=fix(randomu(systime,vals,numbs)*(vals))
          		centbs=fltarr(vals,numbs,2)
              mn=fltarr(numbs)
          		for i=0,vals-1 do centbs[i,0:numbs-1,0:1]=ph[g[i,0:numbs-1],0:1]
               
              ;Standard mean 
              ;meanvar=(moment(mean(centbs[0:vals-1,0:numbs-1,0],dim=1)) )[0:1]

              ;Weighted mean
              wmn=fltarr(numbs)
              FOR i=0,numbs-1 DO BEGIN
                  meanerr2,centbs[0:vals-1,i,0],centbs[0:vals-1,i,1],mx,sx
                  wmn[i]=mx
              ENDFOR
              meanvar=(moment(wmn))[0:1]

          		phase_speed=1./meanvar[0]*delx/delt 
          		speed_error=sqrt(meanvar[1])/meanvar[0]^2*delx/delt

	     ENDIF ELSE BEGIN
             	phase_speed=0
	           	speed_error=0
       ENDELSE
    ENDIF ELSE BEGIN
       	phase_speed=0
   	    speed_error=0
    ENDELSE 


IF n_elements(vals) GT 0 THEN no_tr=vals ELSE no_tr=0

return,[phase_speed,speed_error]

END
