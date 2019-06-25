PRO plot_velvsdw

dates=['20120327','20130914','20130708','20140510','20150122']
bad_ang=list([33*5+5*indgen(7)],[0+5*indgen(5),350,355],[0,350,355],[330,335,340],[0,185+5*indgen(8),345,350,355])

ndates=n_elements(dates)

vel_list=list()
dw_list=list()
FOR i=0,ndates-1 DO BEGIN
	inpath  = 'analysis/CoMP/wave_tracking_output/'+dates[i]+'/'
	restore,inpath+'aa_vel_results.idlsav'
	restore,inpath+'aa_dw_results.idlsav'
	IF n_elements(bad_ang[i]) NE 0 THEN vals=setdifference(fix(results.phi),bad_ang[i],position=good) $
	ELSE good=indgen(n_elements(results))
	vel_list.add,results[good]
	dw_list.add,results_dw[good]  
	
ENDFOR

FOR i=0,ndates-1 DO IF n_elements(vels) EQ 0 THEN vels=reform(vel_list[i].fit_pg[3,*]) ELSE vels=[vels,reform(vel_list[i].fit_pg[3,*])]
FOR i=0,ndates-1 DO IF n_elements(gwids) EQ 0 THEN gwids=reform(vel_list[i].fit_pg[5,*]) ELSE gwids=[gwids,reform(vel_list[i].fit_pg[5,*])]
FOR i=0,ndates-1 DO IF n_elements(dws) EQ 0 THEN dws=dw_list[i].med ELSE dws=[dws,dw_list[i].med]
FOR i=0,ndates-1 DO IF n_elements(delaic) EQ 0 THEN delaic=vel_list[i].aic_line-vel_list[i].aic_pg $ 
			ELSE delaic=[delaic,vel_list[i].aic_line-vel_list[i].aic_pg ]
FOR i=0,ndates-1 DO FOR j=0,n_elements(vel_list[i])-1 DO $
          IF n_elements(pows) EQ 0 THEN pows=total(exp(mypowgauss(alog((vel_list[i].freq)[1:-2,j]),(vel_list[i].fit_pg)[*,j] ))) $
          ELSE pows=[pows,total(exp(mypowgauss(alog((vel_list[i].freq)[1:-2,j]),(vel_list[i].fit_pg)[*,j] ))) ]
FOR i=0,ndates-1 DO FOR j=0,n_elements(vel_list[i])-1 DO $
		  IF n_elements(dat) EQ 0 THEN dat=i ELSE dat=[dat,i]
FOR i=0,ndates-1 DO FOR j=0,n_elements(vel_list[i])-1 DO BEGIN

	p=(vel_list[i].fit_pg)[*,j]
	x=alog((vel_list[i].freq)[1:-2,j])
	par2=double(x-p[4])
	temp=total( p[3]*exp(-(par2/p[5])^2/2.) )
          
    IF n_elements(pows_ga) EQ 0 THEN pows_ga=temp ELSE pows_ga=[pows_ga,temp]
ENDFOR

vels=vels[where(delaic gt 0)]
dws=dws[where(delaic gt 0)]
pows=pows[where(delaic gt 0)]
gwids=gwids[where(delaic gt 0)]
pows_ga=pows_ga[where(delaic gt 0)]
dat=dat[where(delaic gt 0)]

vels=vels[where(dws gt 20)]
pows=pows[where(dws gt 20)]
gwids=gwids[where(dws gt 20)]
pows_ga=pows_ga[where(dws gt 20)]
dat=dat[where(dws gt 20)]
dws=dws[where(dws gt 20)]

dws=sqrt(dws^2-2.*21.^2)

;dws=dws[where(vels gt 1e-4)]
;vels=vels[where(vels gt 1e-4)]
dat_by=bytscl(dat)


fig=scatterplot(dws,vels,xtitle='Non-thermal widths (km s!e-1!n)',ytitle='Enhanced power (km!e2!n s!e-2!n)',$
	            mag=dat_by,rgb_table=13,/ylog,sym_filled=1,symbol='tu',/buffer)
loadct,13,/silent
tvlct,/get,r,g,b
new_ct=transpose([[r],[g],[b]])
text=text(0.75,0.35,vel_list[0,0].date,target=fig)
text=text(0.01,0.88,'b',target=fig,font_size=16,font_style=1)

text=text(0.75,0.31,vel_list[1,0].date,target=fig,color=new_ct[*,dat_by((where(dat eq 1))[0])])
text=text(0.75,0.27,vel_list[2,0].date,target=fig,color=new_ct[*,dat_by((where(dat eq 2))[0])])
text=text(0.75,0.23,vel_list[3,0].date,target=fig,color=new_ct[*,dat_by((where(dat eq 3))[0])])
text=text(0.75,0.19,vel_list[4,0].date,target=fig,color=new_ct[*,dat_by((where(dat eq 4))[0])])
fig.save,'analysis/CoMP/papers/dopp_vel_2/NTW_enhanced_pow.png'

fig2=scatterplot(dws,pows,xtitle='Non-thermal widths (km s!e-1!n)',ytitle='Total power (km!e2!n s!e-2!n)',$
	            mag=dat_by,rgb_table=13,/ylog,sym_filled=1,symbol='tu',/buffer)
loadct,13,/silent
tvlct,/get,r,g,b
new_ct=transpose([[r],[g],[b]])
text=text(0.75,0.35,vel_list[0,0].date,target=fig2)

text=text(0.75,0.31,vel_list[1,0].date,target=fig2,color=new_ct[*,dat_by((where(dat eq 1))[0])])
text=text(0.75,0.27,vel_list[2,0].date,target=fig2,color=new_ct[*,dat_by((where(dat eq 2))[0])])
text=text(0.75,0.23,vel_list[3,0].date,target=fig2,color=new_ct[*,dat_by((where(dat eq 3))[0])])
text=text(0.75,0.19,vel_list[4,0].date,target=fig2,color=new_ct[*,dat_by((where(dat eq 4))[0])])
fig2.save,'analysis/CoMP/papers/dopp_vel_2/NTW_Tot_pow.png'

;Calculate correlation
samps=fix(randomu(systime,n_elements(vels),1000)*n_elements(vels) )
corr_coff=fltarr(1000,2)
corr_coff2=fltarr(1000,2)
FOR i=0,999 DO BEGIN
	corr_coff[i,*]=(r_correlate(dws[samps[*,i]],vels[samps[*,i]]))
	corr_coff2[i,*]=(r_correlate(dws[samps[*,i]],pows[samps[*,i]]))
ENDFOR
print,'Enhanced',mean(corr_coff[*,0]),sqrt((moment(corr_coff[*,0]))[1])
print,'RMS',mean(corr_coff2[*,0]),sqrt((moment(corr_coff2[*,0]))[1])
kstwo,corr_coff,corr_coff2,d,prob
print,'KS',d,prob
stop
END