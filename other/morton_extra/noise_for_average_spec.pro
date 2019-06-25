
FUNCTION noise_for_average_spec,ferr=ferr,er_ch=er_ch,er_ar=er_ar,er_qs=er_qs,see=see

;Temporary variables
date='20120327'
dt=30.

inpath  = 'analysis/CoMP/wave_tracking_output/'+date+'/'

print,n_elements(see)
IF n_elements(see) GT 0 THEN restore,inpath+'err_boxes_seeinc2.sav'$
ELSE restore,inpath+'err_boxes.sav'

sz=size(err)
nx=sz(1)
ny=sz(2)
nt=sz(3)

; define temporal apodization
apodt = fltarr(nt)+1
apod=0.1
apodrimt = nt*apod
apodt[0] = (sin(!pi/2.*findgen(apodrimt)/apodrimt))^2
apodt = apodt*shift(rotate(apodt,2),1) 
CPG=total(apodt)/n_elements(apodt) ;Coherent Power Gain - correction factor needed for power after apodisation

mean_ch=mean(err,dimension=3)
mean_ar=mean(err_ar,dimension=3)
mean_qs=mean(err_qs,dimension=3)




ferr=findgen(nt/2)/(nt*dt)

gnoise=randomn(systime_seed,nx,ny,nt)

noise=fltarr(nx,ny,nt)

FOR i=0,nt-1 DO noise[*,*,i]=gnoise[*,*,i]*mean_ch

err_tot_spec=power_fit_cube(nx,ny,nt,noise,apodt,cpg,noise=1)


res=linfit(findgen(nt/2-1),exp(err_tot_spec[1:nt/2-1,0]))
er_ch=res[0]+res[1]*findgen(nt/2-2)

;AR calculations
sz=size(err_ar)
nx=sz(1)
ny=sz(2)
nt=sz(3)

gnoise=randomn(systime_seed,nx,ny,nt)

noise=fltarr(nx,ny,nt)

FOR i=0,nt-1 DO noise[*,*,i]=gnoise[*,*,i]*mean_ar


err_tot_spec=power_fit_cube(nx,ny,nt,noise,apodt,cpg,noise=1)
res=linfit(findgen(nt/2-1),exp(err_tot_spec[1:nt/2-1,0]))
er_ar=res[0]+res[1]*findgen(nt/2-2)
;window,0
;IF n_elements(SEE) EQ 1 THEN c=1 ELSE plot,err_tot_spec[1:nt/2-1,0]

;QS calculations
sz=size(err_qs)
nx=sz(1)
ny=sz(2)
nt=sz(3)

gnoise=randomn(systime_seed,nx,ny,nt)

noise=fltarr(nx,ny,nt)

FOR i=0,nt-1 DO noise[*,*,i]=gnoise[*,*,i]*mean_qs

err_tot_spec=power_fit_cube(nx,ny,nt,noise,apodt,cpg,noise=1)
res=linfit(findgen(nt/2-1),exp(err_tot_spec[1:nt/2-1,0]))
er_qs=res[0]+res[1]*findgen(nt/2-2)


return,!null


END
