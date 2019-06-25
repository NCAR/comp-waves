;PURPOSE - Takes a data cube and calculates the FFT power for each pixel. 
;          The events with power above 95% significance are kept 
;
;INPUTS - date - date of cube to process
;         dt - cadence of data, if not entered the default is 30s
;         
;OPTIONAL INPUTS sup_pix - creates a super pixel of size given, should probably be used 
;                          for large data cubes for initial test of results
;               set - use to specify which CoMP observable to use, default is velocity
;                     use set='dw' or 'di' to use Doppler width or Intensity
;             
;
;OUTPUTS - probab - provides probability of oscillatory signal is present in each pixel
;                    is probability is above significance. Returns data cube of size n_elements(freq)


PRO comp_waves_fft,date,sup_pix=sup_pix,dt=dt,set=set,probab=probab,freq=freq,$
                   siglvl=siglvl,mx_power=mx_power,cutoff=cutoff,pad=pad,debug=debug
                      


; paths where the CoMP data is located and
; where the .sav-files should be written to
inpath  = 'data/CoMP/wave_tracking_output/'+date+'/'


inpath_ind  = 'data/CoMP/'+strmid(date,8,4,/reverse)+'/'+strmid(date,3,2,/reverse)+'/'+strmid(date,1,2,/reverse)+'/'
files=find_files('*.fts.gz',inpath_ind)
mreadfits,files(0),index

restore,inpath+'cube_ivw_'+date+'.sav' 



;Default is 30 s - may be different for different COMP data sets
IF n_elements(dt) GT 0 THEN dt=dt ELSE dt=30.
IF n_elements(siglvl) EQ 0 THEN siglvl=0.95

IF n_elements(set) EQ 0 THEN set ='dv'

CASE SET OF
   'dv':input=cube_v
   'dw':input=cube_w
   'di':input=cube_i
ENDCASE



sz=size(input)
nx=sz(1)
ny=sz(2)
nt=sz(3)

IF n_elements(sup_pix) EQ 0 THEN sup_pix=1.

nx_r= nx mod sup_pix
ny_r= ny mod sup_pix

dnx=(nx-nx_r)/sup_pix
dny=(ny-ny_r)/sup_pix

Print,'Size of array to be processed'
print,dnx,dny
out=fltarr(dnx,dny)

;set up of variables
rdata=rebin(input(0:nx-nx_r-1,0:ny-ny_r-1,*),dnx,dny,nt)
rdata_mask=rebin(index.mask(0:nx-nx_r-1,0:ny-ny_r-1,*),dnx,dny)

;For padded time-series
IF keyword_set(pad) THEN BEGIN
   IF n_elements(tbas) EQ 0 THEN tbas=3
   base2 = FIX( ALOG(nt) / ALOG(2) ) + tbas   ; power of 2 nearest to N
   padlen=2L^(base2)-nt
   nt=nt+padlen
ENDIF 

freq=findgen(nt/2)/(nt*dt)
probab=fltarr(dnx,dny,nt/2-1)
mx_power=fltarr(dnx,dny,nt/2-1)


;--------------------------------------
;Defines a cut-off value for fitting of power spectra
;Spectra shows break at 0.01 mHz
IF keyword_set(cutoff) THEN cut=nt/2-(where(freq gt 0.01))[0] ELSE cut=0


; define temporal apodization
apodt = fltarr(nt)+1
apod=0.1
apodrimt = nt*apod
apodt[0] = (sin(!pi/2.*findgen(apodrimt)/apodrimt))^2
apodt = apodt*shift(rotate(apodt,2),1) 
CPG=total(apodt)/n_elements(apodt) ;Coherent Power Gain - correction factor needed for power after apodisation

apocube=transpose(rebin(apodt,nt,dny,dnx))

rdata=(rdata-rebin(mean(rdata,dim=3),dnx,dny,nt))*apocube

;DO APODIZATION  
;FOR j=0,dny-1 DO BEGIN
;   FOR i=0,dnx-1 DO BEGIN
;       dum=reform(rdata[i,j,*]) 
;       no_zero=n_elements(where(dum eq 0.))
;       IF no_zero lt 5 THEN BEGIN 
;            rec=dum-mean(dum) 
;            rdata[i,j,*]=rec*apodt
;       ENDIF
;   ENDFOR
;ENDFOR


FOR j=0,dny-1 DO BEGIN
   FOR i=0,dnx-1 DO BEGIN

   dum=reform(rdata[i,j,*])

   if rdata.mask[i,j] EQ 1 THEN BEGIN
      ;Padding of time-series
      s_leno=n_elements(dum)
      IF keyword_set(pad) THEN dum=[dum,fltarr(padlen)]
       
      s_len=n_elements(dum)

      ;APPLY FFT AND CALCULATE POWER
 
     fac=(1.*s_len/s_leno)^2 ; Correction factor required for padding
     power=(2.*abs(fft(dum))^2*fac/CPG^2)[1:nt/2-1-cut]
     
     
     ;Normalised Power
     npw=s_len*power/(moment(dum))[1]
     result = poly_fit(ALOG10(freq(1:nt/2-1-cut)), ALOG10(npw),1) 
    
     IF keyword_set(debug) THEN BEGIN
     	plot,alog10(freq(1:nt/2-1-cut)),alog10(npw)
     	oplot,alog10(freq(1:nt/2-1-cut)),alog10(freq(1:nt/2-1-cut))*result[1]+result[0]
     	pause
     ENDIF

     ; de-bias the normalising factor (see Vaughan 2005)
     N_0 = result[0] + 0.57721466/ALOG(10.0)
     alpha = result[1]
     mod_per = 10.0^N_0 * freq(1:nt/2-1-cut)^alpha
     rat=2.*npw/mod_per
    
     ;Assumes model is correct
     prob = CHISQR_PDF(rat, 2)
     nprob = MAKE_ARRAY(SIZE(prob, /DIM))
     log_pp =  ALOG(DOUBLE(prob));*(n_elements(freq)-1)
     nprob=exp(log_pp)
     
     
      a=where(nprob gt siglvl) 
    

     IF n_elements(a) GT 1 THEN BEGIN
         probab[i,j,a]=nprob[a]
         mx_power[i,j,a]=power[a]
     ENDIF ELSE probab(i,j,*)=0.   

       

   ENDIF
   
   counter,j*dnx+i,dnx*dny,'Calculating',/percent
   ENDFOR
ENDFOR



END