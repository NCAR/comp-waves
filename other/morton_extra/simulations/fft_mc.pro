PRO fft_mc

outpath='analysis/CoMP/simulations/'
IF n_elements(filename) EQ 0 THEN filename='sim_res'

p=[8.5e-6,-1.2,0.001,0.0034,-5.57,0.32]


IF n_elements(fac) EQ 0 THEN fac=1
nt=160*fac ; Length of time-series
nf=nt/2 ; Frequency elements
dt=30./fac   ;Typical cadence


freq = (findgen(nf)+1)/nt/dt
pow=p[0]*freq^p[1];+p[2] ;+p[3]*exp(-(alog(freq)-p[4])^2/2./p[5]^2)
;plot,freq,alog(pow),/xlog,xr=[freq(1),freq(79)]
nelm=150 ; Number of series to averaged together
reps=5000 ;Number of simulation reps
;xarr=fltarr(nt,nelm)
num_bs=500 ;Bootstrapping for SD error on mean
;meanval=fltarr(num_bs)
;varval=fltarr(num_bs)
logfreq=alog(freq[0:nf-2])

local_siglvl=0.95

; define temporal apodization
apodt = fltarr(nt)+1
apod=0.1
apodrimt = nt*apod
apodt[0] = (sin(!pi/2.*findgen(apodrimt)/apodrimt))^2
apodt = apodt*shift(rotate(apodt,2),1) 

mx_freq=fltarr(reps,2)
coeff_fits=fltarr(reps,2)

FOR j=0,reps-1 DO BEGIN


	x = TS_GEN(nt, dt=dt, freq=freq, pow=pow,time=time,/quiet)


	pow_coeff=mle_fit_psd((x-mean(x))*apodt,dt,freqo,pergm)
    coeff_fits[j,*]=pow_coeff

	back_fit_freq=double(10.^(pow_coeff[0]+pow_coeff[1]*ALOG10(freqo)))
    back_fit_freq=back_fit_freq/dt
    pergm_cor=pergm/dt

    ;local_signif=-alog(1-local_siglvl)*back_fit_freq

    ;Signifincae level for multiple measurements
    ;nfc=n_elements(freqo)
	;local_siglvl_for=1-(local_siglvl)^(1./nfc) ;correction for multiple obs
	;global_signif_for =-alog(local_siglvl_for)*back_fit_freq 

    ;mxpow=max(pergm_cor/back_fit_freq,loc)
    ;mxp=freq(loc)
    ;IF pergm_cor[loc] GT local_signif(loc) THEN mx_freq[j,0]=mxp ELSE mx_freq[j,0]=0.
    ;IF pergm_cor[loc] GT global_signif_for(loc) THEN mx_freq[j,1]=mxp ELSE mx_freq[j,1]=0.

    ;plot,freq,pergm_cor,/xlog,/ylog
	;oplot,freq,back_fit_freq
    ;pause

    counter,j,reps,'Calculating',/percent

ENDFOR
stop
END