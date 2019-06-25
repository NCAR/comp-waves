;PURPOSE - Perform simulations to test reliability of CoMP coherence calculations.
;This routine currently looks at the coherence values of two uncorrelated red noise series
;
;Results normally distributed once del_aic lt 0 results removed
;
;

PRO mcsim_coh,just_plot=just_plot,filename=filename,fac=fac,nelm=nelm,coh=coh

outpath='analysis/CoMP/simulations/'
IF n_elements(filename) EQ 0 THEN filename='sim_res_rat'

;FOR plotting
color=['red','blue','green','cyan','black','orange red','dark blue','purple']

IF NOT keyword_set(just_plot) THEN BEGIN
;Average values from Morton 2017
p=[8.5e-6,-1.2,0.001,0.0034,-5.57,0.32]
start=[3e-6,-1.3,5e-3,0.001,-5.5,0.4]

IF n_elements(fac) EQ 0 THEN fac=1
nt=160*fac ; Length of time-series
nf=nt/2 ; Frequency elements
dt=30./fac   ;Typical cadence
freq = (findgen(nf)+1)/nt/dt
pow=p[0]*freq^p[1]+p[2];+p[3]*exp(-(alog(freq)-p[4])^2/2./p[5]^2)
pow2=pow*exp(-100.*freq)
;plot,freq,alog(pow),/xlog,xr=[freq(1),freq(79)]


IF n_elements(nelm) EQ 0 THEN nelm=1 ; Number of series to averaged together
reps=5000 ;Number of simulation reps
xarr=fltarr(nt,nelm)
xarr2=fltarr(nt,nelm)
num_bs=500 ;Bootstrapping for SD error on mean
meanval=fltarr(num_bs)
varval=fltarr(num_bs)

temp={freq:fltarr(nf),coher:fltarr(1),fit:fltarr(2)} 
temp2={spec_ot:complexarr(nf), spec_oo:fltarr(nf), spec_tt:fltarr(nf)}
;Set up filter
fwidth = 0.0015
freq_filt  = 0.0035
nsmooth=8
nspec=nt/2
freq   = findgen(nspec)/(float(nspec*2)*dt)
filter = exp( -(freq-freq_filt)^2/fwidth^2 )
filter(0) = 0.                      ; set dc to zero
filter(where (freq lt .001)) = 0.   ; set low frequencies to zero
filter=filter/total(filter)


;Gain parameter to set series coherence ()
;This formula is derived for white noise - unsure if it holds for red
IF keyword_set(coh) THEN gain=(1.-sqrt(1.-coh^2))/coh ELSE gain=0.0

FOR j=0,reps-1 DO BEGIN

	FOR i=0,nelm-1 DO BEGIN
	    see=2L*i+(nelm)*2L*j ; for selecting different seeds for each time-series
	    
		x = TS_GEN(nt, dt=dt, freq=freq, pow=pow,time=time,/quiet,seed=see)
		xarr[*,i]=x
		see=2L*i+(nelm)*2L*j

		x = TS_GEN(nt, dt=dt, freq=freq, pow=pow,time=time,/quiet,seed=see+1)
		xarr2[*,i]=x
	ENDFOR

	;For creating coherent data
	xarr_temp=xarr+gain*xarr2
	xarr2=xarr2+gain*xarr
	xarr=xarr_temp

    IF nelm EQ 1 THEN BEGIN
    	FT_ser1=fft(xarr-mean(xarr,dim=1),dim=1)
		FT_ser2=fft(xarr2-mean(xarr2,dim=1),dim=1)

    ENDIF ELSE BEGIN
		calc_pow=mean(abs(nt*fft(xarr-transpose(rebin(mean(xarr,dim=1),nelm,nt)),dim=1))^2/2./(nf/2./dt),dim=2)
		calc_pow2=mean(abs(nt*fft(xarr2-transpose(rebin(mean(xarr2,dim=1),nelm,nt)),dim=1))^2/2./(nf/2./dt),dim=2)
	ENDELSE
	
	Xspec_12=smooth(FT_ser1*conj(FT_ser2),nsmooth)
	Xspec_11=real_part(smooth(FT_ser1*conj(FT_ser1),nsmooth))
	Xspec_22=real_part(smooth(FT_ser2*conj(FT_ser2),nsmooth))

	temp.coher=total(filter*abs(Xspec_12)/sqrt(Xspec_11*Xspec_22))
	temp.freq=freq
    temp2.spec_tt=Xspec_22[0:nf-1]
    temp2.spec_oo=Xspec_11[0:nf-1]
    temp2.spec_ot=Xspec_12[0:nf-1]
    ;par=mle_fit_psd_ratio(temp.ratio,30.,freq,nu=2.*nelm)
	;temp.fit=par

	IF n_elements(res) EQ 0 THEN res=temp ELSE res=[temporary(res),temp]
	IF n_elements(res2) EQ 0 THEN res2=temp2 ELSE res2=[temporary(res2),temp2]
				
	counter,j,reps,'How much have we done?',/percent
	
	IF j EQ 100 THEN stop
ENDFOR

stop
save,res,p,pow,nt,dt,nelm,filename=outpath+filename+'.idlsave'
ENDIF

stop

IF keyword_set(just_plot) THEN restore,outpath+filename+'.idlsave'



END	