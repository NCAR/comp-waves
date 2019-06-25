;PURPOSE - Perform simulations to test power law fitting.
;
;Results normally distributed once del_aic lt 0 results removed
;
;

PRO mcsimpowlaw,just_plot=just_plot,filename=filename,fac=fac

outpath='analysis/CoMP/simulations/'
IF n_elements(filename) EQ 0 THEN filename='sim_res'

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
pow=p[0]*freq^p[1]+p[2]+p[3]*exp(-(alog(freq)-p[4])^2/2./p[5]^2)
;plot,freq,alog(pow),/xlog,xr=[freq(1),freq(79)]
nelm=150 ; Number of series to averaged together
reps=5000 ;Number of simulation reps
xarr=fltarr(nt,nelm)
num_bs=500 ;Bootstrapping for SD error on mean
meanval=fltarr(num_bs)
varval=fltarr(num_bs)
logfreq=alog(freq[0:nf-2])
temp={freq:fltarr(nf),spectra:fltarr(nf+1,3),fit_pg:fltarr(6),err_pg:fltarr(6),$
	  chisq_pg:0.0,loglik_pg:0.0,aic_pg:0.0,fit_line:fltarr(3),$
      err_line:fltarr(3),chisq_line:0.0,loglik_line:0.0,aic_line:0.0} 


;Constraints for fitting
;Gaussian model
pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},6)
pi[0].limited(0)=1 & pi[0].limits(0)=1e-8
pi[1].limited(1)=1 & pi[1].limits(1)=-0.1
pi[2].limited(0)=1 & pi[2].limits(0)=5e-6
pi[3].limited(0)=1 & pi[3].limits(0)=1e-5
pi[4].limited([0,1])=1 & pi[4].limits(0)=-7. & pi[4].limits(1)=-5 
pi[5].limited([0,1])=1 & pi[5].limits(0)=-1. & pi[5].limits(1)=1

;Linear model
pi2 = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},3)
pi2[0].limited(0)=1 & pi2[0].limits(0)=1e-8
pi2[1].limited(1)=1 & pi2[1].limits(1)=-0.1
pi2[2].limited(0)=1 & pi2[2].limits(0)=5e-6

FOR j=0,reps-1 DO BEGIN

	FOR i=0,nelm-1 DO BEGIN
		x = TS_GEN(nt, dt=dt, freq=freq, pow=pow,time=time,/quiet)
		xarr[*,i]=x
	ENDFOR


	input=abs(nt*fft(xarr-transpose(rebin(mean(xarr,dim=1),nelm,nt)),dim=1))^2/2./(nf/2./dt)
	indx=round(randomu(systime,nelm,num_bs)*(nelm-1)) 
	
	FOR ix=1,nt/2 DO BEGIN

		logdat=alog(reform(input[ix,*]))
		FOR bi=0,num_bs-1 DO BEGIN
			meanval[bi]=mean(logdat[indx[*,bi]])
			varval[bi]=(moment(logdat[indx[*,bi]]))[1]
		ENDFOR
		temp.spectra[ix,0]=mean(logdat) + 0.57721466 ;Bias correction (Vaughan 2005)
		temp.spectra[ix,1]=sqrt((moment(meanval))[1])
		;ksone,(meanval-temp.spectra[ix,0])/temp.spectra(ix,1),'gauss_cdf',d,prob
		;        		temp.spectra[ix,2]=d 

		IF j EQ 10 THEN BEGIN ;Plots only for rep 10
			IF ix EQ 1 THEN BEGIN
		        plothist,meanval,xbin,pdf,/noplot
		        denestim = akde(meanval,xbin)
		   	    phisto = PLOT(xbin, denestim, XTITLE='Log power', YTITLE='Frequency', $
		           				COLOR=color[0],/buffer,xran=[-3,-8])

		        plothist,logdat,xbin,pdf,/noplot
		        phisto_full = PLOT(xbin, pdf, XTITLE='Log power', YTITLE='Frequency', $
		                		COLOR=color[0],/buffer,/stairstep,xran=[0,-15])

		        plothist,logdat,xbin,pdf,/noplot
		        denestim = akde(logdat,xbin)
		        phisto_full_kde = PLOT(xbin, denestim, XTITLE='Log power', YTITLE='Frequency', $
		                		COLOR=color[0],/buffer,xran=[0,-15])
		                        
		    ENDIF
		                 
		    IF ix MOD 11 EQ 0 THEN BEGIN
		        plothist,meanval,xbin,pdf,/noplot
		        denestim = akde(meanval,xbin)
		    	phisto2 = PLOT(xbin, denestim, XTITLE='Log power', YTITLE='Frequency', $
		                        COLOR=color[ix/10 mod 8],overplot=phisto)

		        plothist,logdat,xbin,pdf,/noplot
		        phisto_full2 = PLOT(xbin, pdf, XTITLE='Log power', YTITLE='Frequency', $
		        COLOR=color[ix/10 mod 8],overplot=phisto_full,/stairstep)

		        plothist,logdat,xbin,pdf,/noplot
		        denestim = akde(logdat,xbin)
		        phisto_full_kde2 = PLOT(xbin, denestim, XTITLE='Log power', YTITLE='Frequency', $
		        COLOR=color[ix/10 mod 8],overplot=phisto_full_kde)
		    ENDIF
		    IF ix EQ nt/2-1 THEN BEGIN
		            phisto.Save, outpath+"/mean_power_plot_MC.jpg"
		            phisto_full.Save, outpath+"/power_plot_MC.jpg"
		            phisto_full_kde.Save, outpath+"/power_plot_MC_kde.jpg"
		    ENDIF   
		ENDIF
	ENDFOR

	respg=mpfitfun('mypowgauss',logfreq,temp.spectra[1:nf-2,0],temp.spectra[1:nf-2,1],start,$
	perror=perror,bestnorm=bestnorm,dof=dof,parinfo=pi,/quiet, nfree=nfree,xtol=1e-5,ftol=1e-5)
	temp.fit_pg=respg
	temp.err_pg=perror
	temp.chisq_pg=bestnorm/dof
	;Calculate log likelihood and Akaike Information Criterion
	temp.loglik_pg=-(nf-2)*alog(2.*!pi)-total(alog(temp.spectra[1:nf-2,1]))-bestnorm/2.
	temp.aic_pg=2.*6.-2*temp.loglik_pg+2.*nfree*(nfree-1)/(dof-1)


	resp=mpfitfun('mypowline',logfreq,temp.spectra[1:nf-2,0],temp.spectra[1:nf-2,1],start[0:2],/quiet,$
	perror=perror2,bestnorm=bestnorm2,dof=dof2,parinfo=pi2,nfree=nfree2,xtol=1e-5,ftol=1e-5)
	temp.fit_line=resp
	temp.err_line=perror2
	temp.chisq_line=bestnorm2/dof2
	;Calculate log likelihood and Akaike Information Criterion
	temp.loglik_line=-(nf-2)*alog(2.*!pi)-total(alog(temp.spectra[1:nf-2,1]))-bestnorm2/2.
	temp.aic_line=2.*nfree2-2*temp.loglik_line+2.*nfree2*(nfree2-1)/(dof2-1)

	IF n_elements(res) EQ 0 THEN res=temp ELSE res=[temporary(res),temp]
				
	counter,j,reps,'How much have we done?',/percent
ENDFOR
save,res,p,pow,nt,dt,nelm,filename=outpath+filename+'.idlsave'
ENDIF

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