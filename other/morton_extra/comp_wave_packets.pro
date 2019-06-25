;PURPOSE - Takes a data cube and calculates the wavelet time-series for each pixel. 
;          The events with power above 95% significance are kept and properties of these
;          significant signals are calculated, e.g. max power, number of wave packets in
;          the signal.
;          Data is de-trended using EMD and largest two IMFs are subtracted as trend
;
;
;INPUTS - date - date of cube to process
;         
;OPTIONAL INPUTS sup_pix - creates a super pixel of size given, should probably be used 
;                          for large data cubes for initial test of results
;               set - use to specify which CoMP observable to use, default is velocity
;                     use set='dw' or 'di' to use Doppler width or Intensity
;               /tc_red - Calculates significance levels for red noise spectra a la Torrence & Compo (1998)
;               /v_red - Calculates significance levels for red noise spectra a la Vaughan (2005)
;               siglvl - significance level for tests - default 0.95
;               
;
;OUTPUTS - mxwv_pack - periods with most wave packets over times series in each pixel
;          out - periods with greatest power in each pixel
;          wv_pack - 3-D array of number of wave packets for each period for each pixel
;


PRO comp_wave_packets,date,sup_pix=sup_pix,local_siglvl=local_siglvl,set=set,$
                      wv_pack=wv_pack,global_siglvl=global_siglvl, $
                      mxwv_pack=mxwv_pack,period=period,out=out,debug=debug
                      


; paths where the CoMP data is located and
; where the .sav-files should be written to
inpath  = 'analysis/CoMP/wave_tracking_output/'+date+'/'
outpath = 'analysis/CoMP/wave_tracking_output/'+date+'/'

restore,inpath+'cube_ivw_'+date+'.sav' 


;Default is 30 s - may be different for different COMP data sets
dt=index.norm_cadence[0]

IF n_elements(siglvl) EQ 0 THEN local_siglvl=0.95
IF n_elements(global_siglvl) EQ 0 THEN global_siglvl=local_siglvl

IF n_elements(set) EQ 0 THEN set ='dv'

CASE SET OF
   'dv':input=cube_v
   'dw':input=cube_w
   'di':input=cube_i
ENDCASE

nx=index.naxis1
ny=index.naxis2
nt=index.nframes

;input=input[*,*,0:nt/2-1]
;nt=(size(input))[3]

;Sets up Super pixel cube
IF n_elements(sup_pix) EQ 0 THEN sup_pix=1.

nx_r= nx mod sup_pix
ny_r= ny mod sup_pix

dnx=(nx-nx_r)/sup_pix
dny=(ny-ny_r)/sup_pix

rdata=rebin(input(0:nx-nx_r-1,0:ny-ny_r-1,*),dnx,dny,nt)
rdata_mask=rebin(index.mask(0:nx-nx_r-1,0:ny-ny_r-1,*),dnx,dny)

Print,'Size of array to be processed'
print,dnx,dny


;-------------------------------------------------------------------
;'Empty' routines to allow set up of variables and confidence intervals

;Find a time-series to test
in=where(rdata_mask eq 1)
test_val=[in[100] mod nx,in[100]/nx]

pow_coeff=mle_fit_psd(reform(rdata[test_val[0],test_val[1],*]),dt,freq,pergm)
nf=n_elements(freq)

wave=wavelet(findgen(nt),dt,period=period,/pad,scale=scale,j=j1,$
                      s0=s0,coi=coi,dj=dj)
;dj=0.1
nscale=n_elements(period)

;Modify to look over certain frequency/period range
p_in=where(period gt 100 and period lt 1000)
f_in=where(freq gt 1./1000 and freq lt 1./100)
coi_in=where(coi gt 100 and coi lt 1000,complement=coi_out)

;Fourier significance
nfc=n_elements(f_in)
local_siglvl_for=1-(global_siglvl)^(1./nfc) ;correction for multiple obs


;Wavelet significance - Auchere et al (2016)
;Corrected for multiple measurements
jcoi=alog(coi[coi_in]/1.033/s0)/alog(2)/dj
nout=total(jcoi>0)
a_coeff = 0.810*(Nout*dj)^0.011
n_coeff = 0.491*(Nout*dj)^0.926
local_siglvl2 = 1 - (1 - global_siglvl^(1/n_coeff))^(1/a_coeff)

Sout = MAX(jcoi)
a_scl_coeff = 0.805+0.45*2^(-Sout*dj)
n_scl_coeff = 1.136*(Sout*dj)^1.2
time_avg_local_siglvl = 1 - (1 - global_siglvl^(1/n_scl_coeff))^(1/a_scl_coeff)
dof = nt - scale/dt

wv_pack=fltarr(dnx,dny,nscale)
mxwv_pack=fltarr(dnx,dny)
mx_per=fltarr(dnx,dny)
mx_freq=fltarr(dnx,dny)

;Define mask for removing COI data
mask=fltarr(nt,nscale,/nozero)
mask[*]=1.
FOR i=0,nt-1 DO BEGIN
   ht=where(period ge coi[i])
   mask(i,ht)=0.
ENDFOR
mask[*,coi_out]=0


;-------------------------------------------------------------------
;Prep for Frequency analysis

; define temporal apodization
apodt = fltarr(nt)+1
apod=0.1
apodrimt = nt*apod
apodt[0] = (sin(!pi/2.*findgen(apodrimt)/apodrimt))^2
apodt = apodt*shift(rotate(apodt,2),1) 
 
;Remove mean  
rdata=(rdata-rebin(mean(rdata,dim=3),dnx,dny,nt))

IF keyword_set(debug) THEN BEGIN
!p.background=255
loadct,5,/silent
window,1
tvim,rdata[*,*,20]<4>(-4)
window,2

ENDIF

;-------------------------------------------------------------------
;Frequency analysis
coeff_fits=fltarr(dnx,dny,2)


FOR j=0,dny-1 DO FOR i=0,dnx-1 DO BEGIN
   scl=0.
   cnt_pix=0.  
   ;Don't calculate for pixels outside mask
   IF rdata_mask[i,j] EQ 1 THEN BEGIN
       rec=reform(rdata[i,j,*]) 
       var=(moment(rec))[1]
             
       wave=wavelet(rec,dt,period=period,/pad,dj=dj,scale=scale,j=j1,$
                      s0=s0,coi=coi)
       ;supresses wave power inside cone of influence - assumes it cannot be trusted
       wpower_orig=abs(wave)^2
       wpower=abs(wave*mask)^2
       power=total(abs(wave*mask)^2,1)/nt ;normalised by variance

       ;Maximum likelihood fitting of power spectra
       ;Following Barret & Vaughan (2012)
       ;base2 = FIX( ALOG(nt) / ALOG(2) ) + 1 
       ;rec_pad = [ rec*apodt, FLTARR( 2L^(base2) - nt ) ]
       pow_coeff=mle_fit_psd(rec*apodt,dt,freq,pergm);,/verbose)
       coeff_fits[i,j,*]=pow_coeff
       ;Calculate power spectra with white noise added
       back_fit_per=double(10.^(pow_coeff[0]+pow_coeff[1]*ALOG10(1./period)))
       back_fit_freq=double(10.^(pow_coeff[0]+pow_coeff[1]*ALOG10(freq)))
       back_fit_per=back_fit_per/dt ;Correct Vaughan normalisation for comparison to wavelet
       back_fit_freq=back_fit_freq/dt
       pergm_cor=pergm/dt
         
          
       ;local significance level
       local_signif=-alog(1-local_siglvl)*back_fit_freq

       ;Following Vaughn (2005) 
       ;Traditional Fourier Method
       global_signif_for =-alog(local_siglvl_for)*back_fit_freq ;Calculation from chi^2 distribution
            
       ;Wavelet global significance
       ;global_signif =-alog(1-local_siglvl2)*back_fit_per
       global_signif=-alog(1-local_siglvl)*back_fit_per
       global_signif = transpose( REBIN(global_signif,j1+1,nt) )

       time_avg_global_signif =-alog(1-time_avg_local_siglvl)*back_fit_per
           
       
          
       ;Calculates frequency with largest power
       ;above background from time averaged wavelet
       ;Only one comparison so use local significance levels
       mxpow=max(power/back_fit_per,loc)
       mxp=period(loc)
       IF power[loc] GT local_signif(loc) THEN mx_per[i,j]=mxp ELSE mx_per[i,j]=0.

       ;Calculates frequency with largest power
       ;above background from time averaged Fourier
       ;Only one comparison so use local significance levels
       mxpow=max(pergm_cor/back_fit_freq,loc)
       mxp=freq(loc)
       IF pergm_cor[loc] GT local_signif(loc) THEN mx_freq[i,j]=mxp ELSE mx_freq[i,j]=0.
    
       ;calculates number of 'wave packets' where power is greater than siglvl significant
       lb=label_region(wpower/global_signif gt 1)
       
       ;go through blobs
      FOR k=1, max(lb) DO BEGIN
          in=where(lb eq k) ; find blob 
          in2=in/nt         ; find which scale each blob row is on
          in3=in mod nt
          pdf=histogram(in2,locations=xbin) ; find number of pixels at each scale
               
          scl=[scl,xbin]
          cnt_pix=[cnt_pix,pdf]
          
      ENDFOR
      
      IF n_elements(scl) GT 1 THEN FOR kk=0,nscale-1 DO wv_pack[i,j,kk]=dt/period[kk]*cnt_pix[where(scl eq kk)]
       ;IF n_elements(a) GT 1 THEN BEGIN
          ;bins number of significant values in each period bin
       ;   res=histogram(a/nt,binsize=1,min=0,max=nscale) 
       ;   wv_pack(i,j,*)=dt/period*res
       ;ENDIF ELSE wv_pack(i,j,*)=0.   

       ;Calculates period with most wave packets
       b=max(wv_pack[i,j,*],loc)
       mxwv_pack[i,j]=period(loc)

       IF keyword_set(debug) THEN BEGIN
          wset,1
          tvim,rdata[*,*,20]<4>(-4)
          plots,i,j,psym=6

          wset,2
          !p.multi=[0,1,2]
          tvim,wpower
          contour,wpower/global_signif,levels=1,/over

          plot,1./freq,pergm_cor,color=0,/xlog,yr=[1e-2,1e3],/ylog,yst=1
          oplot,period,back_fit_per,line=2,color=50
          oplot,1./freq,global_signif_for,line=2,color=10
           
          oplot,period,total(wpower_orig,1)/nt/0.776/0.849,thick=2,color=0
          oplot,period,global_signif[0,*],color=100
          oplot,period,time_avg_global_signif[*],color=150
          xyouts,10.1,8,'Max val'+strtrim(max(wpower/global_signif),2),/data,color=0
                                  

          !p.multi=0
          ;!p.background=0
          ;sdbwave2,rec,delt=30,/fast,/ylog
          ;!p.background=255
          
       ENDIF
        
       counter,j*dnx+i,dnx*dny,'Calculating wavelets',/percent
    ENDIF
ENDFOR


stop
save,coeff_fits,wv_pack,period,filename=outpath+'wave_packets.sav'
stop

loadct,13,/silent
window,0
tvim,out
cgcolorbar,range=[min(out),max(out)]
 
window,1
tvim,mxwv_pack
cgcolorbar,range=[min(mxwv_pack),max(mxwv_pack)]
stop
END