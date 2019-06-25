;PURPOSE - Perform simulations to test power law ratio fitting.
;
;Results normally distributed once del_aic lt 0 results removed
;
;
;Inputs - Fac controls length of time-series, factors of 160
;         nelm controls number of series to average together

PRO mcsimpowlaw_rat,just_plot=just_plot,filename=filename,fac=fac,nelm=nelm

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
dt=30.;/fac   ;Typical cadence

freq = (findgen(nf)+1)/nt/dt
pow=p[0]*freq^p[1]+p[2];+p[3]*exp(-(alog(freq)-p[4])^2/2./p[5]^2)
pow2=pow*exp(-100.*freq)
;plot,freq,alog(pow),/xlog,xr=[freq(1),freq(79)]
IF n_elements(nelm) EQ 0 THEN nelm=1 ; Number of series to averaged together
reps=1;5000 ;Number of simulation reps
xarr=fltarr(nt,nelm)
xarr2=fltarr(nt,nelm)
num_bs=500 ;Bootstrapping for SD error on mean
meanval=fltarr(num_bs)
varval=fltarr(num_bs)

temp={freq:fltarr(nf),ratio:fltarr(nf),fit:fltarr(2),cov:fltarr(2,2)} 



FOR j=0,reps-1 DO BEGIN

	FOR i=0,nelm-1 DO BEGIN
	    see=2L*i+(nelm)*2L*j ; for selecting different seeds for each time-series
	    
		x = TS_GEN(nt, dt=dt, freq=freq, pow=pow,time=time,/quiet,seed=see)
		xarr[*,i]=x
		see=2L*i+(nelm)*2L*j

		x = TS_GEN(nt, dt=dt, freq=freq, pow=pow2,time=time,/quiet,seed=see+1)
		xarr2[*,i]=x
	ENDFOR

    IF nelm EQ 1 THEN BEGIN
    	calc_pow=abs(nt*fft(xarr-mean(xarr,dim=1),dim=1))^2/2./(nf/2./dt)
		calc_pow2=abs(nt*fft(xarr2-mean(xarr2,dim=1),dim=1))^2/2./(nf/2./dt)

    ENDIF ELSE BEGIN
		calc_pow=mean(abs(nt*fft(xarr-transpose(rebin(mean(xarr,dim=1),nelm,nt)),dim=1))^2/2./(nf/2./dt),dim=2)
		calc_pow2=mean(abs(nt*fft(xarr2-transpose(rebin(mean(xarr2,dim=1),nelm,nt)),dim=1))^2/2./(nf/2./dt),dim=2)
	ENDELSE
	
	
	temp.ratio=calc_pow[1:nt/2]/calc_pow2[1:nt/2]
	temp.freq=freq

    par=mle_fit_psd_ratio(temp.ratio,30.,freq,nu=2.*nelm,/plt_lk,cov=cov,filename='/Users/richardmorton/analysis/COMP/simulations/likelihood_test')
	temp.fit=par
	temp.cov=cov

	IF n_elements(res) EQ 0 THEN res=temp ELSE res=[temporary(res),temp]
				
	counter,j,reps,'How much have we done?',/percent
    print, j
	
	
ENDFOR

save,res,p,pow,nt,dt,nelm,filename=outpath+filename+'.idlsave'
ENDIF

stop

IF keyword_set(just_plot) THEN restore,outpath+filename+'.idlsave'

delta_aic=res.aic_line-res.aic_pg
in=where(delta_aic gt 0)
stop
;Histograms of values
pdf=histogram(delta_aic,loc=xbins,binsize=5)
hist=plot(xbins,pdf,/stairstep,/buffer,XTITLE='delta AIC',YTITLE='No. Events')
hist.Save,outpath+'sim_deltaaic_dist'+strtrim(nt,2)+'.png' 

pdf=histogram(alog(res[in].fit_pg[0]),loc=xbins,binsize=0.1)
hist=plot(xbins,pdf,/stairstep,/buffer,XTITLE='Constant',YTITLE='No. Events')
line=plot(alog(p[0])*[1,1],[0,max(pdf)],'--2',overplot=hist)
line=plot(mean(alog(res[in].fit_pg[0]) )*[1,1],[0,max(pdf)],'--2r',overplot=hist)
gf=gaussfit(xbins,pdf,coeff,nterms=3)
line=plot(xbins,gf,'2b',overplot=hist)
hist.Save,outpath+'sim_const_dist'+strtrim(nt,2)+'.png' 

pdf=histogram(res[in].fit_pg[1],loc=xbins,binsize=9e-3)
hist=plot(xbins,pdf,/stairstep,/buffer,XTITLE='Index',YTITLE='No. Events')
line=plot(p[1]*[1,1],[0,max(pdf)],'--2',overplot=hist)
line=plot(mean(res[in].fit_pg[1])*[1,1],[0,max(pdf)],'--2r',overplot=hist)
gf=gaussfit(xbins,pdf,coeff,nterms=3)
line=plot(xbins,gf,'2b',overplot=hist)
hist.Save,outpath+'sim_index_dist'+strtrim(nt,2)+'.png' 

pdf=histogram(res[in].fit_pg[2],loc=xbins, binsize=3e-5)
hist=plot(xbins,pdf,/stairstep,/buffer,XTITLE='Noise',YTITLE='No. Events')
line=plot(p[2]*[1,1],[0,max(pdf)],'--2',overplot=hist)
line=plot(mean(res[in].fit_pg[2])*[1,1],[0,max(pdf)],'--2r',overplot=hist)
gf=gaussfit(xbins,pdf,coeff,nterms=3)
line=plot(xbins,gf,'2b',overplot=hist)
hist.Save,outpath+'sim_noise_dist'+strtrim(nt,2)+'.png' 

pdf=histogram(res[in].fit_pg[3],loc=xbins,binsize=2e-4)
hist=plot(xbins,pdf,/stairstep,/buffer,XTITLE='Amplitude',YTITLE='No. Events')
line=plot(p[3]*[1,1],[0,max(pdf)],'--2',overplot=hist)
line=plot(mean(res[in].fit_pg[3])*[1,1],[0,max(pdf)],'--2r',overplot=hist)
gf=gaussfit(xbins,pdf,coeff,nterms=3)
line=plot(xbins,gf,'2b',overplot=hist)
hist.Save,outpath+'sim_amp_dist'+strtrim(nt,2)+'.png' 

pdf=histogram(res[in].fit_pg[4],loc=xbins,binsize=0.03)
hist=plot(xbins,pdf,/stairstep,/buffer,XTITLE='Central Frequency',YTITLE='No. Events')
line=plot(p[4]*[1,1],[0,max(pdf)],'--2',overplot=hist)
line=plot(mean(res[in].fit_pg[4])*[1,1],[0,max(pdf)],'--2r',overplot=hist)
gf=gaussfit(xbins,pdf,coeff,nterms=3)
line=plot(xbins,gf,'2b',overplot=hist)
hist.Save,outpath+'sim_cent_freq_dist'+strtrim(nt,2)+'.png'

pdf=histogram(abs(res[in].fit_pg[5]),loc=xbins,binsize=0.02)
hist=plot(xbins,pdf,/stairstep,/buffer,XTITLE='Width',YTITLE='No. Events')
line=plot(p[5]*[1,1],[0,max(pdf)],'--2',overplot=hist)
line=plot(mean(abs(res[in].fit_pg[5]) )*[1,1],[0,max(pdf)],'--2r',overplot=hist)
gf=gaussfit(xbins,pdf,coeff,nterms=3)
line=plot(xbins,gf,'2b',overplot=hist)
hist.Save,outpath+'sim_width_dist'+strtrim(nt,2)+'.png'  

END	