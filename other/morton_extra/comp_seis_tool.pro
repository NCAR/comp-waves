
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
;HISTORY: Created by RJM JUNE 2016
;
;
;

FUNCTION comp_seis_tool, vel, delx, delt, debug=debug, ret=ret,no_tr=no_tr,algn_err=algn_err,$
                          gauss=gauss,vel_orig=vel_orig,ph=ph,errors=errors



COMMON located_dat,located, nx, nt
COMMON threads_dat, threads

s=size(vel)
ntime=s(1)
nzs=s(2)


;errors=sqrt(mean( (vel-smooth(vel,[3,1],/edge_truncate) )^2,dim=1)) ; Estimated error on velocity
IF n_elements(vel_orig) LT 1 THEN vel_orig=vel
;atrous,vel,decomp=a,n_scales=1
vel=smooth(vel,[3,1],/edge_truncate);a[*,*,0] ;+smooth(a[*,*,1],3,/edge_truncate)
;vel=a[*,*,0]

seis_res=create_struct('s0',fltarr(1,2))

IF n_elements(algn_err) GT 0 THEN print,'Error on alignment is:', algn_err


numg=0.001 ;gradient cutoff
numt=5     ; minimum thread length
meas_size=9 


    IF not keyword_set(gauss) THEN BEGIN
      ;output is in common located_dat
      locate_things_min,data=vel[*,*],grad=numg
      
      ;print,'no Gaussian fit used'
    ENDIF ELSE BEGIN
      IF n_elements(errors) gt 0. THEN errorsi=transpose(rebin(errors,nzs,ntime))
      
      locate_things_fg,data=vel[*,*],grad=numg,meas_size=meas_size, $
                 errors=errorsi,check=check,cut_chisq=cut_chisq
    
      IF n_elements(algn_err) GT 0 THEN located.errs=located.errs+algn_err   ; adds alignment error to error from Gaussian fits

      ;print,'Gaussian fit used'
    ENDELSE

    window,0
    !p.multi=[0,2,1]
    tvim,vel
    tvim,located.peaks[*,*,1]
    
    ;output is in common thread_dat
    follow_thread_fg,min_tlen=numt,area=3
    threads_max=threads

    ;Do for minimum strips
    IF not keyword_set(gauss) THEN BEGIN
      ;output is in common located_dat
      locate_things_min,data=vel[*,*]*(-1.),grad=numg
      
      ;print,'no Gaussian fit used'
    ENDIF ELSE BEGIN
      IF n_elements(errors) gt 0. THEN errorsi=transpose(rebin(errors,nzs,ntime))
      
      locate_things_fg,data=vel[*,*]*(-1.),grad=numg,meas_size=meas_size, $
                 errors=errorsi,check=check,cut_chisq=cut_chisq
    
      IF n_elements(algn_err) GT 0 THEN located.errs=located.errs+algn_err   ; adds alignment error to error from Gaussian fits

      ;print,'Gaussian fit used'
    ENDELSE
    
    ;output is in common thread_dat
    follow_thread_fg,min_tlen=numt,area=3
    threads_mn=threads
   
    
    IF n_elements(threads_mn) gt 1 THEN threads=[threads_max,threads_mn[1:n_elements(threads_mn)-1]] $
    ELSE threads=threads_max

    
    sz=size(threads)
    len=0
    IF sz[1] GT 1 THEN BEGIN
     	 ph=fltarr(sz[1]-1,15)
       h=1
     	 FOR i=1,sz[1]-1 DO BEGIN
        	dummy=threads[i].pos
        	terr=threads[i].err_pos
          dummy_vel=threads[i].inten
          terr_vel=threads[i].err_inten
        	;if value missing set to same as last pixel and set
        	zer=where(dummy EQ 0.)
        	IF zer[0] NE -1 THEN BEGIN
           		FOR ii=0,n_elements(zer)-1 DO BEGIN
                    		dummy[zer[ii]]=dummy[zer[ii]-1]
                    		terr[zer[ii]]=1.
                        dummy_vel[zer[ii]]=dummy_vel[zer[ii]-1]
                        terr_vel[zer[ii]]=1.
           		ENDFOR
        	ENDIF
        	in=where(dummy gt -1.)
          IF float(n_elements(zer))/n_elements(in) LT 0.3 THEN BEGIN
              	;len=[len,n_elements(in)]
                ;seis_res=create_struct(temporary(seis_res),'s'+strtrim(h,2), $
                 ;                      [[in],[dummy[in]],[interpolate(vel,dummy[in],in,cubic=-0.5)],[errors[in]]] )
                 seis_res=create_struct(temporary(seis_res),'s'+strtrim(h,2), $
                                       [[in],[dummy[in]],[dummy_vel[in]],[terr_vel[in]]] )
                h=h+1 ;counter
                IF keyword_set(debug) THEN BEGIN
                   !p.multi=[0,2,2]
                   plot,in,interpolate(vel,dummy[in],in,cubic=-0.5)
                   oploterr,in,interpolate(vel,dummy[in],in,cubic=-0.5),errors[in]
                    
                   plot,in,dummy[in],yst=1
                   grad=deriv(dummy[in]*delt,in*delx)
                  print,grad
                   oplot,in,min(dummy[in])+mean(grad)*findgen(n_elements(in)),linestyle=2

                   plot,in,grad
                   plot,in,/nodata
                   pause
                ENDIF
                !p.multi=[0,2,1]
                loadct,13,/silent
                window,1
                plot,findgen(n_elements(in)),dummy[in],xst=1,yst=1
                oploterr,findgen(n_elements(in)),dummy[in],terr[in]

              	res=poly_fit(findgen(n_elements(in)),dummy[in],1,yfit=fit,measure_err=terr[in],sigma=sig,chisq=chi)
              	ph[i-1,0]=res[1] ;gradient
                ph[i-1,1]=sig[1] ;error
                ph[i-1,2]=chi/(n_elements(in)-2) ;reduced chi^2
                ph[i-1,3]=n_elements(in)         ;length of track
                ph[i-1,4]=in[0]
                oplot,findgen(n_elements(in)),fit,linestyle=1,color=100
                res=poly_fit(findgen(n_elements(in)),dummy[in],2,yfit=fit,measure_err=terr[in],sigma=sig,chisq=chi)
                ph[i-1,5]=res[1]
                ph[i-1,6]=res[2]
                ph[i-1,7]=sig[1] ;error
                ph[i-1,8]=sig[2] ;error
                ph[i-1,9]=chi/(n_elements(in)-3) ;reduced chi^2
                oplot,findgen(n_elements(in)),fit,linestyle=2,color=200
                res=mpfitfun('myexp',findgen(n_elements(in)),dummy[in],terr[in],[min(dummy[in]),0.1],perror=sig,bst=chi,/quiet)
                ph[i-1,10]=res[0]
                ph[i-1,11]=res[1]
                ph[i-1,12]=sig[0] ;error
                ph[i-1,13]=sig[1] ;error
                ph[i-1,14]=chi
                oplot,findgen(n_elements(in)),res[0]*exp(findgen(n_elements(in))*res[1]),linestyle=3
                print,reform(ph[i-1,*])
                
                plot,dummy_vel[in]
                oploterr,dummy_vel[in],terr_vel[in]
                pause
            ENDIF
	     ENDFOR
       !p.multi=0
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
