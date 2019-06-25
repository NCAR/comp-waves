;PURPOSE - Perform simulations to test reliability of CoMP coherence calculations.
;This time looking at how coherence behaves with Welch method

;This routine currently looks at the coherence values of two uncorrelated red noise series
;
;
;

PRO mc_sim_coh_pt2,just_plot=just_plot,filename=filename,fac=fac,nelm=nelm,coh=coh

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
pow=p[0]*freq^p[1]+p[2]+0.01;+p[3]*exp(-(alog(freq)-p[4])^2/2./p[5]^2)

;plot,freq,alog(pow),/xlog,xr=[freq(1),freq(79)]


IF n_elements(nelm) EQ 0 THEN nelm=1 ; Number of series to averaged together
reps=5000 ;Number of simulation reps
tser=fltarr(nt,nelm)
tser2=fltarr(nt,nelm)

;Welch parameters
length=60
overlap=0.2
ps_welch,findgen(nt),length=length,overlap=overlap,ft_ser=ft_ser1 ;test

sz=size(ft_ser1)
nspec=sz(1)/2

temp={freq:fltarr(nspec),spec_ot:complexarr(nspec,sz(2)), spec_oo:fltarr(nspec,sz(2)), spec_tt:fltarr(nspec,sz(2)), $
	  tseries:fltarr(nt,nelm), tseries2:fltarr(nt,nelm)}

freq_welch   = findgen(nspec)/(float(nspec*2)*dt)



;Gain parameter to set series coherence ()
;This formula is derived for white noise - unsure if it holds for red
IF keyword_set(coh) THEN gain=(1.-sqrt(1.-coh^2))/coh ELSE gain=0.0

FOR j=0,reps-1 DO BEGIN

	FOR i=0,nelm-1 DO BEGIN
	    see=2L*i+(nelm)*2L*j ; for selecting different seeds for each time-series
	    
		;x = TS_GEN(nt, dt=dt, freq=freq, pow=pow,time=time,/quiet,seed=see)
		tser[*,i]=randomn(see,nt) ;x
		see=2L*i+(nelm)*2L*j

		;x = TS_GEN(nt, dt=dt, freq=freq, pow=pow,time=time,/quiet,seed=see+1)
		tser2[*,i]=randomn(see+1,nt)  ;x
	ENDFOR
    temp.tseries=tser
    temp.tseries2=tser2

	;For creating coherent data
	tser_temp=tser+gain*tser2
	tser2=tser2+gain*tser
	tser=tser_temp

    ps_welch,tser,length=length,overlap=overlap,ft_ser=ft_ser1
    ps_welch,tser2,length=length,overlap=overlap,ft_ser=ft_ser2

	Xspec_12=FT_ser1*conj(FT_ser2)
	Xspec_11=real_part(FT_ser1*conj(FT_ser1))
	Xspec_22=real_part(FT_ser2*conj(FT_ser2))

	
	temp.freq=freq_welch
    temp.spec_tt=Xspec_22[0:nspec-1,*]
    temp.spec_oo=Xspec_11[0:nspec-1,*]
    temp.spec_ot=Xspec_12[0:nspec-1,*]
    ;par=mle_fit_psd_ratio(temp.ratio,30.,freq,nu=2.*nelm)
	;temp.fit=par

	IF n_elements(res) EQ 0 THEN res=temp ELSE res=[temporary(res),temp]
	
				
	counter,j,reps,'How much have we done?',/percent
	
	;IF j EQ 100 THEN stop
ENDFOR

stop
save,res,p,pow,nt,dt,nelm,filename=outpath+filename+'.idlsave'
ENDIF

stop

IF keyword_set(just_plot) THEN restore,outpath+filename+'.idlsave'



END	