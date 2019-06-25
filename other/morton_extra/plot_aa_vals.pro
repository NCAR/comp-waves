; For comparison of results from aa_comp_vals.pro
;
; One task is to perform MC simulations of fitting procedure
; to look for parameter correlations
;
;
;INPUTS date - vector of dates to restore
;		ang_remove - list containing arrays of angles 
;                    of regions that are considered bad (see aa_comp_vel angle plots)
;                    If no regions need removing add empty array to list
;
;e.g., IDL>dates=['20120327','20130914','20130708','20140510','20150121']
;      IDL> bad_ang=list([33*5+5*indgen(7)],[0+5*indgen(5),350,355],[0,350,355],[330,335,340],[0,185+5*indgen(8),345,350,355])
;
;
;CALLS: setdifference from CG library
;
;

;Freedman and Diaconis bin width estimate
FUNCTION calc_bin_wdith,data

     nelm=n_elements(data)
 	 r=iqr(data)
 	 w=2.*r/nelm^(1./3.)

 	 return,w
END	

FUNCTION create_table,date_dict,date,outpath

	openw,lun2,outpath+'fit_table.txt',/get_lun
	printf,lun2,'Date & Angle & A & n & B & C & D & E & $\chi_{\nu}^2$\\'
	printf,lun2,'     &       & ($10^{-6}$) &   & ($10^{-3}$) & ($10^{-3}$) &  &  & \\'
	printf,lun2,'\hline \\'
    FOR i=0,n_elements(date)-1 DO BEGIN
    	datep=date_dict[date[i]].date
    	angle=date_dict[date[i]].phi
    	vals=date_dict[date[i]].fit_pg
    	err=date_dict[date[i]].err_pg
    	chi=date_dict[date[i]].chisq_pg  
    	delta_AIC=date_dict[date[i]].aic_line-date_dict[date[i]].aic_pg;+6*4./(71.-3.-1.)-12*7./(71.-6.-1.)
    	FOR j=0,n_elements(date_dict[date[i]])-1 DO BEGIN
    	IF delta_AIC[j] GT 0 THEN $
    	printf,lun2,datep[0]+' & '+strtrim(string(angle[j],format='(-F10.0)'),2) $
    	+' & '+strtrim(string(vals[0,j]*1e6,format='(-F10.2)'),2)+'$\pm$'+strtrim(string(err[0,j]*1e6,format='(-F10.2)'),2) $
    	+' & '+strtrim(string(vals[1,j],format='(-F10.2)'),2)+'$\pm$'+strtrim(string(err[1,j],format='(-F10.2)'),2) $
    	+' & '+strtrim(string(vals[2,j]*1e3,format='(-F10.2)'),2)+'$\pm$'+strtrim(string(err[2,j]*1e3,format='(-F10.2)'),2) $
    	+' & '+strtrim(string(vals[3,j]*1e3,format='(-F10.2)'),2)+'$\pm$'+strtrim(string(err[3,j]*1e3,format='(-F10.2)'),2) $
    	+' & '+strtrim(string(vals[4,j],format='(-F10.2)'),2)+'$\pm$'+strtrim(string(err[4,j],format='(-F10.2)'),2) $
    	+' & '+strtrim(string(abs(vals[5,j]),format='(-F10.2)'),2)+'$\pm$'+strtrim(string(err[5,j],format='(-F10.2)'),2) $
    	+' & '+strtrim(string(chi[j],format='(-F10.2)'),2)+'\\'
    	ENDFOR
    ENDFOR
    printf,lun2,'\hline \\'
	close,lun2
	free_lun,lun2
END


PRO plot_aa_vals,date,ang_remove=ang_remove,domc=domc,print2file=print2file

resolve_routine, 'gauss_cdf',/is_function

ndates=n_elements(date)
outpath='analysis/CoMP/papers/dopp_vel_2/'
IF keyword_set(print2file) THEN openw,lun,outpath+'key_vaules.txt',/get_lun

;Read in data
date_dict=hash()
FOR i=0,ndates-1 DO BEGIN
	inpath  = 'analysis/CoMP/wave_tracking_output/'+date[i]+'/'
	restore,inpath+'aa_vel_results.idlsav'

	;Manually remove poor regions
	IF n_elements(ang_remove[i]) NE 0 THEN vals=setdifference(fix(results.phi),ang_remove[i],position=good) $
	ELSE good=indgen(n_elements(results))

	p=median(results[good].fit_pg,dim=2)
	IF keyword_set(print2file) THEN printf,lun,'Median values  '+date[i],p
    date_dict[date[i]]=results[good]
ENDFOR


;FOR plotting
color=['red','blue','green','cyan','black','orange red','dark blue','purple']

;freq
freq=list((date_dict[date[0]])[0].freq)
FOR i=1,ndates-1 DO freq.add,(date_dict[date[i]])[0].freq

;AIC model comparison
delta_AIC=fltarr(72,ndates)
;Added correction for n/k<=0 - should add to main code.
FOR i=0,ndates-1 DO delta_AIC[0,i]=date_dict[date[i]].aic_line-date_dict[date[i]].aic_pg;+6*4./(71.-3.-1.)-12*7./(71.-6.-1.)
pdf=histogram(delta_AIC[where(delta_aic NE 0)],loc=xbins,binsize=10)
aic_hist=plot(xbins,pdf,/stairstep,/buffer,XTITLE='$\Delta_{AIC}$',YTITLE='No. Events',xran=[-20,300])
aic_hist.Save,outpath+'aic_hist.png'


;slope
slope_fits=fltarr(72,ndates)
FOR i=0,ndates-1 DO slope_fits[0,i]=date_dict[date[i]].fit_pg[1]

in=where(slope_fits NE 0 and delta_AIC gt 0)

pdf=list()
xbins=list()
pdf.add,histogram(slope_fits[in],loc=xbin,binsize=0.05)
xbins.add,xbin
FOR i=0,ndates-1 DO BEGIN
	in2=where(slope_fits[*,i] NE 0 and delta_AIC[*,i] gt 0)
	pdf.add,histogram(slope_fits[in2,i],loc=xbin,binsize=0.1)
	xbins.add,xbin
	print,n_elements(slope_fits[in2,i])
ENDFOR

mn=moment(slope_fits[in])
med=median(slope_fits[in])
IF keyword_set(print2file) THEN printf,lun,'Slope fits',mn,med,n_elements(slope_fits[in])


name=list()
for i=0,ndates-1 DO name.add,strmid(date[i],6,2)+'-'+strmid(date[i],4,2)+'-'+strmid(date[i],0,4)

slope_hist=plot(xbins[0],pdf[0],/stairstep,/buffer,XTITLE='Power law index',YTITLE='No. Events',name='All',thick=2)
slope_hist1=plot(xbins[1],pdf[1],'r',/stairstep,overplot=slope_hist,transparency=40,thick=2,name=name[0])
slope_hist2=plot(xbins[2],pdf[2],'royal blue',/stairstep,overplot=slope_hist,transparency=40,thick=2,name=name[1])
slope_hist3=plot(xbins[3],pdf[3],'lime green',/stairstep,overplot=slope_hist,transparency=40,thick=2,name=name[2])
slope_hist4=plot(xbins[4],pdf[4],'dark orange',/stairstep,overplot=slope_hist,transparency=40,thick=2,name=name[3])
slope_hist5=plot(xbins[5],pdf[5],'deep pink',/stairstep,overplot=slope_hist,transparency=40,thick=2,name=name[4])
line=plot(mn[0]*[1,1],[0,max(pdf[0])],'--2r',overplot=slope_hist)
line=plot(med[0]*[1,1],[0,max(pdf[0])],'--2b',overplot=slope_hist)
leg=legend(target=[slope_hist,slope_hist1,slope_hist3,slope_hist2,slope_hist4,slope_hist5],position=[-0.1,29],/data,$
			auto_text=1)
text=text(0.02,0.9,'d',font_size=16, font_style=1,font_name='Helvetica')
slope_hist.Save,outpath+'slope_hist.png'

;Scale factor for power spectra
const_prop=fltarr(72,ndates)
FOR i=0,ndates-1 DO const_prop[0,i]=date_dict[date[i]].fit_pg[0]

IF keyword_set(print2file) THEN printf,lun,'const',moment(const_prop[in]),median(const_prop[in])



;maximum frequency
max_freq=fltarr(72,ndates)
FOR i=0,ndates-1 DO BEGIN 
	
	len=n_elements(date_dict[date[i]])
	FOR j=0,len-1 DO BEGIN
		part1=(date_dict[date[i]])[j].fit_pg[0]*freq[i]^(date_dict[date[i]])[j].fit_pg[1]+(date_dict[date[i]])[j].fit_pg[2]

		z=(alog(freq[i])-(date_dict[date[i]])[j].fit_pg[4])^2/2./(date_dict[date[i]])[j].fit_pg[5]^2
		part2=(date_dict[date[i]])[j].fit_pg[3]*exp(-z)
		mx=max(part2/part1,loc)
		max_freq[j,i]=freq[i,loc]
		
	ENDFOR
ENDFOR

mn=moment(max_freq[in])
med=median(max_freq[in])
IF keyword_set(print2file) THEN printf,lun,'Max Frequency',mn,med
pdf=histogram(max_freq[in],loc=xbins,binsize=5e-4)
maxfreq_hist=plot(xbins,pdf,/stairstep,/buffer,XTITLE='Maximum frequency (Hz)',YTITLE='No. Events',xran=[0,0.01])
line=plot(mn[0]*[1,1],[0,max(pdf)],'--2r',overplot=maxfreq_hist)
line=plot(med[0]*[1,1],[0,max(pdf)],'--2b',overplot=maxfreq_hist)
maxfreq_hist.Save,outpath+'maxfreq_hist.png'




;central frequency
cent_freq=fltarr(72,ndates)
FOR i=0,ndates-1 DO cent_freq[0,i]=date_dict[date[i]].fit_pg[4]
mn=moment(exp(Cent_freq[in]))
med=median(exp(Cent_freq[in]))
IF keyword_set(print2file) THEN printf,lun,'Cent Frequency',mn,med

pdf=list()
xbins=list()
pdf.add,histogram(exp(cent_freq[in]),loc=xbin,binsize=5e-4)
xbins.add,xbin
FOR i=0,ndates-1 DO BEGIN
	in2=where(cent_freq[*,i] NE 0 and delta_AIC[*,i] gt 0)
	pdf.add,histogram(exp(cent_freq[in2,i]),loc=xbin,binsize=5e-4)
	xbins.add,xbin
ENDFOR

;Run hypothesis test
in_freq=where(cent_freq[*,2] NE 0 and delta_AIC[*,2] gt 0)
vals1=exp(cent_freq[in_freq,2])
FOR i=0,ndates-1 DO BEGIN
	 in_freq=where(cent_freq[*,i] NE 0 and delta_AIC[*,i] gt 0)
     vals2=exp(cent_freq[in_freq,i])
     IF i eq 2 THEN IF keyword_set(print2file) THEN printf,lun,'nout' ELSE IF keyword_set(print2file) THEN printf,lun,'Wilcoxon',2.*(rs_test(vals1,vals2))[1] ;two-tailed
     kstwo,vals1,vals2,D,prob
     IF keyword_set(print2file) THEN printf,lun,'K-S 2',prob
     
ENDFOR

cent_hist=plot(xbins[0]*1e3,pdf[0],/stairstep,/buffer,XTITLE='Central frequency (mHz)',$
			  YTITLE='No. Events',xran=[0,0.01]*1e3,name='All', thick=2)
cent_hist1=plot(xbins[1]*1e3,pdf[1],'r',/stairstep,overplot=cent_hist,transparency=40,thick=2,name=name[0])
cent_hist2=plot(xbins[2]*1e3,pdf[2],'royal blue',/stairstep,overplot=cent_hist,transparency=40,thick=2,name=name[1])
cent_hist3=plot(xbins[3]*1e3,pdf[3],'lime green',/stairstep,overplot=cent_hist,transparency=40,thick=2,name=name[2])
cent_hist4=plot(xbins[4]*1e3,pdf[4],'dark orange',/stairstep,overplot=cent_hist,transparency=40,thick=2,name=name[3])
cent_hist5=plot(xbins[5]*1e3,pdf[5],'deep pink',/stairstep,overplot=cent_hist,transparency=40,thick=2,name=name[4])
;plotg=plot(freq(0),80*exp(-(freq(0)-mn[0])^2/2./mn[1]),overplot=cent_hist)
leg=legend(target=[cent_hist,cent_hist1,cent_hist3,cent_hist2,cent_hist4,cent_hist5],position=[0.0095*1e3,75],/data,$
			auto_text=1)
line=plot((mn[0])*[1,1]*1e3,[0,max(pdf[0])],'--2r',overplot=cent_hist)
line=plot((med[0])*[1,1]*1e3,[0,max(pdf[0])],'--2b',overplot=cent_hist)
text=text(0.02,0.9,'a',font_size=16, font_style=1,font_name='Helvetica')
cent_hist.Save,outpath+'cent_hist.png'

;Gaussian width
cent_width=fltarr(72,ndates)
FOR i=0,ndates-1 DO cent_width[0,i]=date_dict[date[i]].fit_pg[5]

mn=moment(abs(cent_width[in]))
med=median(abs(cent_width[in]))
IF keyword_set(print2file) THEN printf,lun,'Width',mn,med


pdf=list()
xbins=list()
pdf.add,histogram(abs(cent_width[in])/alog(10),loc=xbin,binsize=0.03)
xbins.add,xbin
FOR i=0,ndates-1 DO BEGIN
	in2=where(abs(cent_width[*,i]) NE 0 and delta_AIC[*,i] gt 0)
	pdf.add,histogram(abs(cent_width[in2,i])/alog(10),loc=xbin,binsize=0.03)
	xbins.add,xbin
ENDFOR

gw_hist=plot(xbins[0],pdf[0],/stairstep,/buffer,XTITLE='Characteristic width in frequency decades',$
			  YTITLE='No. Events',xran=[0,0.5],name='All', thick=2)
gw_hist1=plot(xbins[1],pdf[1],'r',/stairstep,overplot=gw_hist,transparency=40,thick=2,name=name[0])
gw_hist2=plot(xbins[2],pdf[2],'royal blue',/stairstep,overplot=gw_hist,transparency=40,thick=2,name=name[1])
gw_hist3=plot(xbins[3],pdf[3],'lime green',/stairstep,overplot=gw_hist,transparency=40,thick=2,name=name[2])
gw_hist4=plot(xbins[4],pdf[4],'dark orange',/stairstep,overplot=gw_hist,transparency=40,thick=2,name=name[3])
gw_hist5=plot(xbins[5],pdf[5],'deep pink',/stairstep,overplot=gw_hist,transparency=40,thick=2,name=name[4])
;plotg=plot(freq(0),80*exp(-(freq(0)-mn[0])^2/2./mn[1]),overplot=gw_hist)
leg=legend(target=[gw_hist,gw_hist1,gw_hist3,gw_hist2,gw_hist4,gw_hist5],position=[0.45,75],/data,$
			auto_text=1)
line=plot((mn[0])*[1,1]/alog(10),[0,max(pdf[0])],'--2r',overplot=gw_hist)
line=plot((med[0])*[1,1]/alog(10),[0,max(pdf[0])],'--2b',overplot=gw_hist)
text=text(0.02,0.9,'b',font_size=16, font_style=1,font_name='Helvetica')

;pdf=histogram(abs(cent_width[in])/alog(10),loc=xbins,binsize=0.02)
;gw_hist=plot(xbins,pdf,/stairstep,/buffer,XTITLE='Peak width in frequency decades',YTITLE='No. Events')
;line=plot(mn[0]/alog(10)*[1,1],[0,max(pdf)],'--2r',overplot=gw_hist)
;line=plot(med[0]/alog(10)*[1,1],[0,max(pdf)],'--2b',overplot=gw_hist)
gw_hist.Save,outpath+'gw_hist.png'


;Peak amplitude
peak_amp=fltarr(72,ndates)
FOR i=0,ndates-1 DO peak_amp[0,i]=date_dict[date[i]].fit_pg[3]/$
exp(mypowline(date_dict[date[i]].fit_pg[4],date_dict[date[i]].fit_pg[0:2]))

mn=moment(abs(peak_amp[in]))
med=median(abs(peak_amp[in]))
IF keyword_set(print2file) THEN printf,lun,'Peak amp',mn,med


pdf=list()
xbins=list()
pdf.add,histogram(alog10(peak_amp[in]),loc=xbin,binsize=0.2)
xbins.add,xbin
FOR i=0,ndates-1 DO BEGIN
	in2=where(abs(peak_amp[*,i]) NE 0 and delta_AIC[*,i] gt 0)
	pdf.add,histogram(alog10(peak_amp[in2,i]),loc=xbin,binsize=0.2)
	xbins.add,xbin
ENDFOR

pamp_hist=plot(xbins[0],pdf[0],/stairstep,/buffer,XTITLE='log!d10!n Relative Peak Power',$
			  YTITLE='No. Events',xran=[-2,2],name='All', thick=2)
pamp_hist1=plot(xbins[1],pdf[1],'r',/stairstep,overplot=pamp_hist,transparency=40,thick=2,name=name[0])
pamp_hist2=plot(xbins[2],pdf[2],'royal blue',/stairstep,overplot=pamp_hist,transparency=40,thick=2,name=name[1])
pamp_hist3=plot(xbins[3],pdf[3],'lime green',/stairstep,overplot=pamp_hist,transparency=40,thick=2,name=name[2])
pamp_hist4=plot(xbins[4],pdf[4],'dark orange',/stairstep,overplot=pamp_hist,transparency=40,thick=2,name=name[3])
pamp_hist5=plot(xbins[5],pdf[5],'deep pink',/stairstep,overplot=pamp_hist,transparency=40,thick=2,name=name[4])

leg=legend(target=[pamp_hist,pamp_hist1,pamp_hist3,pamp_hist2,pamp_hist4,pamp_hist5],position=[-0.2,55],/data,$
			auto_text=1)
line=plot(alog10(mn[0])*[1,1],[0,max(pdf[0])],'--2r',overplot=pamp_hist)
line=plot(alog10(med[0])*[1,1],[0,max(pdf[0])],'--2b',overplot=pamp_hist)
text=text(0.02,0.9,'c',font_size=16, font_style=1,font_name='Helvetica')


pamp_hist.Save,outpath+'pamp_hist.png'

;Effective SD
eff_sd=fltarr(100,69,ndates)
FOR i=0,ndates-1 DO eff_sd[0,0,i]=date_dict[date[i]].spectra[1:-1,1]
in2=where(eff_sd ne 0)
pdf=histogram(eff_sd[in2],loc=xbins,binsize=0.01)
aic_hist=plot(xbins,pdf,/stairstep,/buffer,XTITLE='Effective $\sigma$ (km!e2!n s!e-2!n)',YTITLE='No. Events')
line=plot(!pi/sqrt(6.*150)*[1,1],[0,2000],'--2',overplot=aic_hist)
aic_hist.Save,outpath+'eff_sd_hist.png'

;Relative excess
pow=list(0)
pow_ll=list(0)
pow_gg=list(0)
;72 hardcoded for 5 degree angle slices in corona
FOR i=0,ndates-1 DO pow.add,fltarr(n_elements(freq[i]),72)
FOR i=0,ndates-1 DO pow_ll.add,fltarr(n_elements(freq[i]),72)
FOR i=0,ndates-1 DO pow_gg.add,fltarr(n_elements(freq[i]),72)
pow_dat=list(0)
FOR i=0,ndates-1 DO pow_dat.add,fltarr(n_elements(freq[i]),72)

noise=fltarr(72,ndates)
FOR i=0,ndates-1 DO noise[0,i]=date_dict[date[i]].fit_pg[2]
IF keyword_set(print2file) THEN printf,lun,'Noise',moment(noise[in]),median(noise[in])
amp=fltarr(72,ndates)
FOR i=0,ndates-1 DO amp[0,i]=date_dict[date[i]].fit_pg[3]
IF keyword_set(print2file) THEN printf,lun,'Max pow',moment(amp[in]),median(amp[in])

excess_pow=fltarr(71,ndates)
FOR j=0,ndates-1 DO BEGIN 
	phi=date_dict[date[j]].phi
	FOR i=0,71 DO BEGIN
	
	IF n_elements(where(phi eq i*5)) EQ 1 THEN $
	IF delta_AIC[i,j] GT 0 THEN BEGIN
	    nf=n_elements(freq[j])
		pow_l=const_prop[i,j]*freq[j,1:nf-1]^slope_fits[i,j]+noise[i,j]
		pow_g=pow_l+amp[i,j]*exp(-((alog(freq[j,1:nf-1]))-cent_freq[i,j])^2/2./cent_width[i,j]^2)
	    pow[j+1,1:nf-1,i]=pow_g/pow_l-1.
        pow_ll[j+1,1:nf-1,i]=pow_l
		pow_gg[j+1,1:nf-1,i]=pow_g
        pow_dat[j+1,1:nf-1,i]=exp((date_dict[date[j]].spectra[1:-1,0])[*,i])/pow_l-1
        excess_pow[i,j]=int_tabulated(freq[j,1:nf-1],pow_g)/int_tabulated(freq[j,1:nf-1],pow_l)
        
    ENDIF
    ENDFOR
    
ENDFOR

totpow=fltarr(72,ndates)
FOR i=0,ndates-1 DO totpow[0,i]=mean(pow[i+1,1:-1,*],dim=1)

pdf=histogram(totpow[where(delta_AIC[0:71,*] gt 0)],loc=xbins,binsize=5)
totpow_hist=plot(xbins,pdf,/stairstep,/buffer,XTITLE='Relative excess power',$
				YTITLE='No. Events',xran=[0,200],xst=1)
totpow_hist.Save,outpath+'tot_pow_hist.png'

av_pow=list(0)
med_pow=list(0)
av_pow_fit=list(0)
FOR i=0,ndates-1 DO BEGIN
    in3=where(totpow[*,i] GT 0)
    mn=moment(pow_dat[i+1,1:-1,in3],dim=2)
    av_pow.add,mn
    med_pow.add,mean(pow_dat[i+1,1:-1,in3],dim=2)
    av_pow_fit.add,moment(pow[i+1,1:-1,in3],dim=2)
ENDFOR

plotex=plot(freq[0,1:-1]*1e3,med_pow[1,1:-1,0],line=6,symbol='plus',/buffer,xtitle='Frequency (mHz)',$
	ytitle='Normalised excess power',color='red',name=name[0])
plotex2a=plot(freq[1,1:-1]*1e3,med_pow[2,1:-1,0],color='royal blue',line=6,symbol='tu',overplot=plotex,name=name[1])
plotex3a=plot(freq[2,1:-1]*1e3,med_pow[3,1:-1,0],color='lime green',line=6,symbol='s',overplot=plotex,name=name[2])
plotex4a=plot(freq[3,1:-1]*1e3,med_pow[4,1:-1,0],color='dark orange',line=6,symbol='p',overplot=plotex,name=name[3])
plotex5a=plot(freq[4,1:-1]*1e3,med_pow[5,1:-1,0],color='deep pink',line=6,symbol='o',overplot=plotex,name=name[4])
plotex1b=plot(freq[0,1:-1]*1e3,av_pow_fit[1,1:-1,0],'2',color='red',overplot=plotex)
plotex2b=plot(freq[1,1:-1]*1e3,av_pow_fit[2,1:-1,0],'2',color='royal blue',overplot=plotex)
plotex3b=plot(freq[2,1:-1]*1e3,av_pow_fit[3,1:-1,0],'2',color='lime green',overplot=plotex)
plotex4b=plot(freq[3,1:-1]*1e3,av_pow_fit[4,1:-1,0],'2',color='dark orange',overplot=plotex)
plotex5b=plot(freq[4,1:-1]*1e3,av_pow_fit[5,1:-1,0],'2',color='deep pink',overplot=plotex)
leg=legend(target=[plotex,plotex2a,plotex3a,plotex4a,plotex5a],position=[0.018*1e3,2.35],/data)
text=text(0.02,0.9,'b',font_size=16, font_style=1,font_name='Helvetica')
plotex.Save,outpath+'rel_diff.png'

res=create_table(date_dict,date,outpath)

;Plot against MDI data
;For PSD correction we have multiplied by s_1^2/f_s s_2 (the random looking number)
files=find_files('*.fts','analysis/mdi_144day/')
mdi_pow=readfits(files)
mdi_freq=findgen(1296)*6.43e-6
mdi_pow=total(mdi_pow,1)
nsmo=1
mdi_comp1=plot(date_dict[date[0]].freq*1e3,exp(smooth(mean(date_dict[date[0]].spectra[0:-2,0],dim=2),nsmo))*30*153.67,/ylog,/xlog,$
	          xr=[1e-4,3e-2]*1e3,yr=[1,6e2],thick=2,xtit="Frequency (mHz)",ytit='Velocity PSD (km!e2!ns!e-2!n Hz!e-1!n)',color='red',$
	          name=name[0],/buffer,yst=1)

mdi_comp2=plot(date_dict[date[1]].freq*1e3,exp(smooth(mean(date_dict[date[1]].spectra[0:-2,0],dim=2),nsmo))*30*142.56,$
	           thick=2,overplot=mdi_comp1,  name=name[1],color='royal blue',line=2)
mdi_comp3=plot(date_dict[date[2]].freq*1e3,exp(smooth(mean(date_dict[date[2]].spectra[0:-2,0],dim=2),nsmo))*30*154*1.133,$
			thick=2,color='lime green',overplot=mdi_comp1, name=name[2],line=2)
mdi_comp4=plot(date_dict[date[3]].freq*1e3,exp(smooth(mean(date_dict[date[3]].spectra[0:-2,0],dim=2),nsmo))*30*168.48,$
			thick=2,'-.',color='dark orange',overplot=mdi_comp1, name=name[3],line=3)
mdi_comp5=plot(date_dict[date[4]].freq*1e3,exp(smooth(mean(date_dict[date[4]].spectra[0:-2,0],dim=2),nsmo))*30*168.48,$
			thick=2,'-',color='Gold',overplot=mdi_comp1, name=name[4],line=4)
;mdi_comp=plot(mdi_freq,(mdi_pow/max(mdi_pow))^0.8,overplot=mdi_comp1,transparency=50,color='blue violet', $
;				name='MDI Doppler')
mdi_comp5.order,/send_to_back
;plot some measure of variability
mn=moment(date_dict[date[1]].spectra[0:-2,0],dim=2)
upper=smooth(exp(mn[*,0]),nsmo)+exp(mn[*,0])*sqrt(mn[*,1]/65)
lower=smooth(exp(mn[*,0]),nsmo)-exp(mn[*,0])*sqrt(mn[*,1]/65)
;line1=errorplot([0,2e-4],[0,3.5],[0,(upper[75]-lower[75])*30*142.56],overplot=mdi_comp1,thick=2,symbol='dot')
;line1=plot([2e-4,2e-4],3+[0,1]*(upper[75]-lower[75])*30*142.56,overplot=mdi_comp1,symbol='tu',thick=2)
;line1=plot(date_dict[date[1]].freq[0:-2],lower*30*142.56,overplot=mdi_comp1)
restore,'analysis/comp/2007/coronal_average_specrebin_1.idlsav'
mdi_comp6=plot(spec2007.freq[4:-2]*1e3,smooth(exp(spec2007.spec[4:-2,0])*29.*324.93,1),overplot=mdi_comp1,$
	line=5,thick=2,name='30-10-2005')
;mdi_comp_mf=plot(exp(mean(Cent_freq[in]))*[1,1],[2e-4,1e-1],transparency=50,thick=2,overplot=mdi_comp)
leg=legend(target=[mdi_comp6,mdi_comp1,mdi_comp3,mdi_comp2,mdi_comp4,mdi_comp5],position=[0.9,0.85])
text=text(0.05,0.88,'a',target=mdi_comp1,font_size=16,font_style=1)
mdi_comp1.Save,outpath+'total_corona_vs_mdi.png'

;##########################################################
;Fit average coronal power spectra
pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},6)
pi[0].limited(0)=1 & pi[0].limits(0)=1e-8
pi[1].limited(1)=1 & pi[1].limits(1)=-0.1
pi[2].limited(0)=1 & pi[2].limits(0)=5e-6
pi[3].limited(0)=1 & pi[3].limits(0)=1e-5
pi[4].limited([0,1])=1 & pi[4].limits(0)=-7. & pi[4].limits(1)=-5 
pi[5].limited(1)=1 & pi[5].limits(1)=0.8

startvals=double([2e-6,-1.2,5e-3,0.001,-5.5,0.4])

corr=[30*153.67,30*142.56,30*154*1.133,30*168.48,30*168.48,29.*947.93]

FOR i=0,ndates-1 DO BEGIN
    spectra2f=mean(date_dict[date[i]].spectra[1:-2,0],dim=2);+corr[i]
    x2f=alog(date_dict[date[i]].freq[1:-2])
    errors2f=sqrt((moment(date_dict[date[i]].spectra[1:-2,0],dim=2))[*,1])
	respg=mpfitfun('mypowgauss',x2f[*,0],spectra2f,errors2f,startvals,$
	perror=perror,bestnorm=bestnorm,dof=dof,parinfo=pi,/quiet,nfree=nfree)
    IF keyword_set(print2file) THEN printf,lun,'Average fit values', respg
    IF keyword_set(print2file) THEN printf,lun,'Average fit values error', perror
    
ENDFOR

spectra2f=spec2007.spec[4:-2,0];+corr[5]
x2f=alog(spec2007.freq[4:-2])
errors2f=spec2007.spec[4:-2,1];+corr[5]
respg=mpfitfun('mypowgauss',x2f[*,0],spectra2f,errors2f,startvals,$
perror=perror,bestnorm=bestnorm,dof=dof,parinfo=pi,/quiet,nfree=nfree)
IF keyword_set(print2file) THEN printf,lun,'Average fit values', respg
IF keyword_set(print2file) THEN printf,lun,'Average fit values error', perror

;##########################################################
;plot comparing AIA & CoMP data
restore,'analysis/sdo/2012/03/27/20120327_1600-2000_north_ch_15Mm_mean_nuwt_vel_amp.sav'
nuwt_va_ch=nuwt_vel_amp
restore,'analysis/sdo/2012/03/27/20120327_1600-2000_east_qs_50Mm_mean_nuwt_vel_amp.sav'
nuwt_va_qs=nuwt_vel_amp
restore,'analysis/comp/wave_tracking_output/20120327/boxed_average_velocity_power_spectra.idlsav'
npfit_ch=read_csv('analysis/r_table_cut.csv',header=thead)
npfit_qs=read_csv('analysis/r_table_20120327_qs_50mm.csv',header=thead_qs)


trend_ch=smooth(exp(tot_spec[1:-1,0])*30*153.67,3)
trend_qs=smooth(exp(mean((date_dict[date[0]])[50:52].spectra[*,0],dim=2))*30*153.67,3,/edge_trunc)
err_qs=trend_qs*mean((date_dict[date[0]])[50:52].spectra[*,1],dim=2)
aiacomp=plot(f(1:-1)*1e3,trend_ch,/xlog,/ylog,xtitle='Frequency (mHz)',$
	ytitle='Velocity PSD (arbitrary)',yr=[1e-1,0.5e3],xr=[9e-5,3e-2]*1e3,yst=1,xst=1,name='CoMP CH',thick=2,/buffer,$
	font_size=16,position=[0.16,0.1,0.92,0.92])
aiacomp_err=plot(f(1:-1)*1e3,trend_ch+tot_spec[1:-1,2]*30*153.67,line=2,overplot=aiacomp)
aiacomp_err=plot(f(1:-1)*1e3,trend_ch-tot_spec[1:-1,2]*30*153.67,line=2,overplot=aiacomp)
aiacomp_qs=plot(date_dict[date[0]].freq*1e3,trend_qs/10.,overplot=aiacomp,color='magenta',name='CoMP QS',thick=2)
aiacomp_err=plot(date_dict[date[0]].freq*1e3,trend_qs/10.+err_qs/10.,line=2,overplot=aiacomp,color='magenta')
aiacomp_err=plot(date_dict[date[0]].freq*1e3,trend_qs/10.-err_qs/10.,line=2,overplot=aiacomp,color='magenta')


erb=npfit_ch.field3-npfit_ch.field4
cis=[exp(npfit_ch.field3+0.577)*erb]
res=poly_fit(alog(npfit_ch.field2[10:37]),npfit_ch.field3[10:37]+0.577,1,measure_err=erb[10:37]+0.577,sigma=sigma)

aiacomp3=errorplot(npfit_ch.field2[0:37]*1e3,exp(npfit_ch.field3[0:37]+0.577)/1.2,cis[0:37]/1.2,'tu',color='dodger blue',overplot=aiacomp,linesty=6,$
					errorbar_color='dodger_blue',name='SDO AIA CH',/sym_filled)
;dr=plot(npfit_ch.field2[10:38],10^(res[0]/alog(10))*npfit_ch.field2[10:38]^(res[1]),overplot=aiacomp,color='red')
aiacomp3.order,/send_to_back

erb2=npfit_qs.field3-npfit_qs.field4
cis2=[exp(npfit_qs.field3+0.577)*erb2]
res=poly_fit(alog(npfit_qs.field2[10:37]),npfit_qs.field3[10:37]+0.577,1,measure_err=erb2[10:37]+0.577,sigma=sigma)

aiacomp4=errorplot(npfit_qs.field2[0:37]*1e3,exp(npfit_qs.field3[0:37]+0.577)/190,cis2[0:37]/190,'s',color='lime green',overplot=aiacomp,linesty=6,$
					errorbar_color='lime green',name='SDO AIA QS',/sym_filled)
aiacomp4.order,/send_to_back
;dr=plot(npfit_qs.field2[10:37],10^(res[0]/alog(10))*npfit_qs.field2[10:37]^(res[1])/170,overplot=aiacomp,color='red')
leg=legend(target=[aiacomp,aiacomp3,aiacomp_qs,aiacomp4],pos=[0.89,0.85])

text=text(0.01,0.88,'d',target=aiacomp,font_size=18,font_style=1)

aiacomp.Save,outpath+'aia_vs_comp.png'


;##########################################################
aiain=where(nuwt_va_ch.peak_freq lt 0.03)
aiaamp=plot(nuwt_va_ch.peak_freq[aiain]*1e3,nuwt_va_ch.peak_vel_amp[aiain],'+',/xlog,/ylog,xtitle='Frequency (mHz)',$
	ytitle='Velocity Amplitude (km s!e-1!n)',linesty=6,position=[0.12,0.1,0.62,0.6],xst=1,yst=1,yran=[1,100],$
	 xran=[1e-4,0.1]*1e3,/buffer,font_size=14)
pdf = HISTOGRAM(nuwt_va_ch.peak_vel_amp[aiain], LOCATIONS=xbin,binsize=2)
b = PLOT(xbin,pdf, FILL_COLOR='red',/current,position=[0.67,0.1,0.92,0.6],xran=[100,1],/xlog,/stairstep, $
		xshowtext=0,fill_background=1,fill_trans=30,font_size=14)
b.rotate,-90
b.scale,1.64,0.6
b.histogram=1
pdf = HISTOGRAM(nuwt_va_ch.peak_freq[aiain], LOCATIONS=xbin,binsize=5e-4)
b2 = PLOT(xbin*1e3,pdf, FILL_COLOR='blue',/current,position=[0.12,0.65,0.62,0.9],xran=[1e-4,0.1]*1e3,/xlog,/stairstep, $
		xshowtext=0,fill_background=1,fill_trans=30,font_size=14)
b2.histogram=1
text=text(0.01,0.88,'c',target=aiaamp,font_size=16,font_style=1)
aiaamp.Save,outpath+'aia_waves.png'


aiain=where(nuwt_va_qs.peak_freq lt 0.03)
aiaamp2=plot(nuwt_va_qs.peak_freq[aiain]*1e3,nuwt_va_qs.peak_vel_amp[aiain],'tu',color='green',$
			 xst=1,yst=1,yran=[1,100],xran=[1e-4,0.1]*1e3,/buffer,/xlog,/ylog,xtitle='Frequency (mHz)',$
	ytitle='Velocity Amplitude (km s!e-1!n)',linesty=6,position=[0.1,0.1,0.6,0.6])
pdf = HISTOGRAM(nuwt_va_qs.peak_vel_amp[aiain], LOCATIONS=xbin,binsize=2)
b = PLOT(xbin,pdf, FILL_COLOR='red',/current,position=[0.65,0.1,0.9,0.6],xran=[100,1],/xlog,/stairstep, $
		xshowtext=0,fill_background=1,fill_trans=30)
b.rotate,-90
b.scale,1.64,0.6
b.histogram=1
pdf = HISTOGRAM(nuwt_va_qs.peak_freq[aiain], LOCATIONS=xbin,binsize=5e-4)
b2 = PLOT(xbin/1e3,pdf, FILL_COLOR='blue',/current,position=[0.1,0.65,0.6,0.9],xran=[1e-4,0.1]/1e3,/xlog,/stairstep, $
		xshowtext=0,fill_background=1,fill_trans=30)
b2.histogram=1
aiaamp2.Save,outpath+'aia_waves2.png'

stop

IF keyword_set(print2file) THEN	close,lun
IF keyword_set(print2file) THEN	free_lun,lun

stop
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
restore,'analysis/CoMP/wave_tracking_output/20120327/flat_limb.sav'
restore,'analysis/CoMP/wave_tracking_output/20120327/cube*'
out=rotate(out,1)
sz=size(out)
nx=sz(1)
ny=sz(2)
;shift by 90 degrees
out=shift(out,90*nx/360)
;rebin
scal_fact=2.8
out=congrid(out,nx/scal_fact,ny)

;get aia colour table
aia_lct,wav='193',/load
tvlct,red,green,blue,/get
cols=[[red],[green],[blue]]

rsun=695.7

flat_plot=image(out,rgb_table=cols,axis_style=0,layout=[1,4,1])
xaxis=axis('x',location=0,coord_trans=[0,360./nx*scal_fact],title='Angle (theta)')
yaxis=axis('y',location=0,coord_trans=[(220)*index.xscale/rsun,index.xscale/rsun],title='Distance (R_{sun})')

grad_plot=plot(date_dict['20120327'].phi,date_dict['20120327'].fit_pg[1],'2d-',xtitle='Angle (theta)',$
				ytitle='Power law Gradient',layout=[1,4,2],/current,xst=1)

freq_plot=plot(date_dict['20120327'].phi,1e3*exp(date_dict['20120327'].fit_pg[4]),'2d-',xtitle='Angle (theta)',$
				ytitle='Frequency (mHz)',layout=[1,4,3],/current,xst=1)

amp_plot=plot(date_dict['20120327'].phi,date_dict['20120327'].fit_pg[3],'2d-',xtitle='Angle (theta)',$
				ytitle='Power' ,layout=[1,4,4],/current,xst=1)



stop
END