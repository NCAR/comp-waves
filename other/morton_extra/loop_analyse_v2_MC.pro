;LOOP ANALYSE.pro
;
;PURPOSE: Used to analyse coronal loops from CoMP data. Will separate inward and outward
;propagating power for both loop legs to enable estimate of damping & propagation speeds.
;
;
;
;OPTIONAL INPUTS - time_split - if set will break the time-series in two to
;                               calculate power spectra for each part of series
;                  foot - two element array specifying footpoints
;                  slice - scalar specifying middle slice - default is middle
;                  n_slice - scalar number of additional slices to use - default is 10
;
;HISTORY - Written by RJM 2015
;
;POTENTIAL THINGS TO DO: * When calculating fits to lags, all values of lag coherence are considered 
;                          Maybe filter out 'bad' lags at least using null hypothesis test.
;

;#########################################################################################
FUNCTION loop_pow_mc,data,loop_num,cph_in,cph_out,lgno,outpath,inpow,outpow,ratio_mc,$
                  freq_n,k_n,halfk

COMMON loop_param,nx,nz,nt,cadence,pix,mid,sc,ns,mid_eo
COMMON fourier_param,tbas,ybas,half,freq

IF lgno EQ 1 THEN rot=4 ELSE rot=1 

;set up array for storing results
TD_diagram=transpose( reform(data[0,0:-1,0:-1]) ) 
power=twod_fft(TD_diagram,cadence,pix,freq_n=freq_n,k_n=k_n,tbas=tbas,ybas=ybas,/han)
sz=size(power)
nyt=sz[1] & ntt=sz[2]
half=(ntt-1)/2
halfk=(nyt-1)/2
em_corr= 0.57721566 ;Euler constant

;Should always be even now
;IF ntt mod 2 eq 0 THEN lbound=half+2 ELSE lbound=half+1
;IF ntt mod 2 eq 0 THEN sub=1 ELSE sub=0

nslice=ns[1]-ns[0]+1
midr=(size(TD_diagram))[1]

twod_pow=fltarr((nyt-1)/2+1,ntt,nslice)

;Loop over slices
FOR i=ns[0],ns[1] DO BEGIN

    TD_diagram=transpose(reform(data[i,0:-1,0:-1])) 
    FOR j=0,midr-1 DO TD_diagram[j,0:-1]=TD_diagram[j,0:-1]-mean(TD_diagram[j,0:-1])
     
    power=twod_fft(TD_diagram,cadence,pix,tbas=tbas,ybas=ybas,/han)
    twod_pow[1:halfk,0:half-1,i-ns[0]]=reverse(power[1:halfk,1:half],2)
    twod_pow[1:halfk,half+1:ntt-2,i-ns[0]]=reverse(power[1:halfk,half+2:ntt-1],2)
    IF n_elements(temp) EQ 0 THEN temp=power ELSE temp=[[[temp]],[[power]]]
ENDFOR
twod_pow=twod_pow[0:-1,0:ntt-2,0:-1]

  
;Just for plotting
ks=reverse(freq_n(0:half+1))/cph_in
ks_out=reverse(freq_n(0:half+1))/cph_out  

;Summed power
;All power below 0.05 Mm-1, i.e. as in Verth et al.,
in=where(k_n lt 0.08)

;IF n_elements(twod_pow[0,*,0]) mod 2 EQ 0 THEN midt=[half+2,2*half+1] ELSE midt=[half+1,2*half]

;Sum over wavenumber, then find mean value over
;the different contributing slices
;Work in log-space, with bias correction (em_corr)
IF rot eq 4 THEN BEGIN
   inpow=reverse(alog(total(twod_pow[1:max(in),0:half-1,*],1))+em_corr)
   outpow=alog(total(twod_pow[1:max(in),half+1:2*half,*],1))+em_corr
ENDIF ELSE BEGIN
   outpow=reverse(alog(total(twod_pow[1:max(in),0:half-1,*],1))+em_corr)
   inpow=alog(total(twod_pow[1:max(in),half+1:2*half,*],1))+em_corr
ENDELSE


;Find mean and SD via bootstrap
num_bs=500
indx=floor(randomu(systime,nslice,num_bs)*(nslice))
ratio=inpow-outpow    ;log-space so subtract
ratio_mean=mean(ratio,dim=2)
varval=fltarr(half,num_bs) 
FOR bi=0,num_bs-1 DO BEGIN
     varval[*,bi]=(mean(ratio[0:-1,indx[0:-1,bi]],dim=2))
ENDFOR

var=sqrt((moment(varval,dim=2))[0:-1,1])
stop


;##############################################
;Plot mean loop power
;Obtain mean values for each of the five hundred BSMC samples
;Appears approximately normal (maybe log-normal)
new=(moment(twod_pow,dim=3,/double))[*,*,0:1]
new=(moment(temporary(new[*,*,*,0]),dim=3,/double))[*,*,0:1]

kw=[[[rotate(new[*,*,0],rot)]],[[rotate(new[*,*,1],rot)]]]

;removing NAN's
kw[*,0,*]=0.
kw[half+1,*,*]=0.
kw[0,*,*]=0.

loadct,13,/silent
device,/encapsul,/color,filename=outpath+'power_'+strtrim(loop_num,2)+'st_half_leg_'+strtrim(lgno,2)+'.eps'

vals=moment(alog10(kw[*,*,0])>(-10)<0.)
bot=-4.7 ;vals[0]-sqrt(vals[1])
top=-2.7 ;vals[0]+sqrt(vals[1])


tvim,alog10(kw[*,*,0])<(top)>(bot),yrange=[0,k_n(halfk)],xrange=[-freq_n[half],freq_n[half]]*1e3,aspect=2.5,$
                          xtitle='Frequency (mHz)',ytitle=' Wavenumber (Mm!e-1!n)'

loadct,0,/silent
oplot,-reverse(freq_n(1:half))*1e3,ks(1:half),color=255,thick=2
oplot,freq_n(1:half)*1e3,reverse(ks_out(1:half)),color=255,thick=2

device,/close

return,kw
END


;#########################################################################################
;#########################################################################################
;Calculates the propagation of velocity signals via cross-correlation
;Uses noise estimated from signal to produce bootstrap repetitions of
;CC for each row.
;Borrowed and modified from compute_speed.pro by S. Tomczyk

FUNCTION calc_speed, n,mid,nt,icent,pro_vel,errs,num_mc,cor=cor

;##############################################################
;Cross-correlate adjacent rows in time-distance map
;##############################################################
cent=fltarr(n)
FOR i=-n/2,n/2 DO BEGIN
  mlag=ccpeak(pro_vel(*,icent),pro_vel(*,icent+i))
  cent(i+n/2)=mlag
ENDFOR

;##############################################################
;Calculate derivative of the lags giving 1/c_ph & take median value
;##############################################################
c=[0.,median(deriv(cent))]   ;use robust estimate for mean shift;  collapse all time series along track
col=fltarr(nt)


;##############################################################
;Align, assuming constant speed, to create an improved S/N time-series
;##############################################################
FOR i=0,mid-1 DO BEGIN
  shift=c(1)*(i-mid/2)+c(0)
  col=col+interpolate(pro_vel(*,i),findgen(nt)+shift)
ENDFOR
col=col/float(mid)


;##############################################################
;Use collapsed time series to perform cross correlation with all time series along track
;##############################################################
cent=fltarr(mid,num_mc)
cor=fltarr(mid,num_mc)
dcent=fltarr(mid)+.1

FOR j=0,num_mc-1 DO BEGIN
    FOR i=0,mid-1 DO BEGIN
        mlag=ccpeak(col,smooth(pro_vel(*,i),3,/edge_truncate)+errs(*,i,j),c_cor=ccor)
        cent(i,j)=mlag
        cor(i,j)=ccor
   ENDFOR
ENDFOR

return,cent
END


;#########################################################################################
;#########################################################################################
;INPUTS - num_mc - number of Monte Carlo data sets to generate
;       - pos- loop positions to work with, e.g., footpoint to apex
;       - n_slice - number of slices to use
;
;OPTIONAL INPUTS - /second_leg - set if using second half of loop
;
;OUTPUTS - lag_pro/ret - lag values for prograde/retrograde waves
;          fits_p/fits_r - results of linear fits to lag_pro/ret
;          cph - array of phase speeds and errors
;

PRO wave_vel_calc,data,num_mc,pos,n_slice,second_leg=second_leg,lags_pro=lags_pro,$
                  lags_ret=lags_ret,fits_p=fits_p,fits_r=fits_r,cph=cph,resid_r=resid_r,resid_p=resid_p

COMMON loop_param,nx,nz,nt,cadence,pix,mid,sc,ns,mid_eo

midr=pos[1]-pos[0]+1

icent=midr/2
nlag=nt  ;number of points in cross correlation (make odd)
lag=indgen(nlag)-fix(nlag/2.)
n = midr + ((midr mod 2)-1)  ;number of points along track for initial guess (make odd)



IF n_slice mod 2 eq 0 THEN cslic=n_slice+1 ELSE cslic=n_slice
randerr=randomn(systime,nt,midr,num_mc) ;generate random errors
lags_pro=fltarr(midr,num_mc,cslic,2)
lags_ret=fltarr(midr,num_mc,cslic,2)

IF keyword_set(second_leg) THEN rot=7 ELSE rot=0


print,'Starting bootstrap calculation - make take a while if num_mc is big'

FOR k=ns[0],ns[1] DO BEGIN

    vel=rotate(reform(data[k,0:nt-1,pos[0]:pos[1]]),rot)
    FOR i=0,midr-1 DO vel[*,i]=vel[*,i]-mean(vel[*,i])

    nmid=(midr-1)/2+1
    tmid=(nt-1)/2+1
    ;odd 0:nmid-1 or nmid:end
    ;even 0:nmid-1 or nmid+1:end
    IF midr mod 2 EQ 0 THEN bound=nmid+1 ELSE bound=nmid
    IF nt mod 2 EQ 0 THEN tbound=tmid+1 ELSE tbound=tmid

    trans=fft(vel)
    pro_trans=trans             ;select prograde waves (assume nt even, npt odd)
    pro_trans(0:tmid-1,0:nmid-1)=0.
    pro_trans(tbound:nt-1,bound:midr-1)=0.
    pro_vel=float(fft(pro_trans,/inverse))
    
    ;estimate noise from data - uses variation at highest spatial frequency
    noise=fltarr(midr)
    FOR i=0,midr-1 DO noise[i]=sqrt(mean((pro_vel[*,i]-smooth(pro_vel[*,i],3,/edge_truncate))^2))
    
    
    noises=randerr*transpose(rebin(noise,midr,nt,num_mc),[1,0,2])
    lags_pro[*,*,k-ns[0],0]=calc_speed(n,midr,nt,icent,pro_vel,noises,num_mc,cor=cor)
    lags_pro[*,*,k-ns[0],1]=cor

  
    ret_trans=trans             ;select retrograde waves
    ret_trans(0:tmid-1,bound:midr-1)=0.
    ret_trans(tbound:nt-1,0:nmid-1)=0.
    ret_vel=float(fft(ret_trans,/inverse))

    FOR i=0,midr-1 DO noise[i]=sqrt(mean((ret_vel[*,i]-smooth(ret_vel[*,i],3,/edge_truncate))^2))
    
    noises=randerr*transpose(rebin(noise,midr,nt,num_mc),[1,0,2])
    lags_ret[*,*,k-ns[0],0]=calc_speed(n,midr,nt,icent,ret_vel,noises,num_mc,cor=cor) 
    lags_ret[*,*,k-ns[0],1]=cor
    
ENDFOR


;set_up arrays to store variables
fits_p=fltarr(num_mc,5)
fits_r=fltarr(num_mc,5)
resid_p=fltarr(midr,num_mc)
resid_r=fltarr(midr,num_mc)

dist=(findgen(midr)-midr/2)*pix


FOR i=0,num_mc-1 DO BEGIN
     lagstf=reform(lags_pro[*,i,*,0])
     corval=reform(lags_pro[*,i,*,1])
     ;good_p=where(corval gt 0.4) ;Where cross-correlation values less than 0.4    

     mn=(moment(lagstf,dimension=2))[*,0:1] 

     ;For some reason, 0 value of mn is always large and positive!
     ;So leave out of fit
     res=poly_fit(dist[1:midr-1],mn[1:midr-1,0]*cadence,1,measure=sqrt(mn[1:midr-1,1])*cadence,sigma=sig,chisq=chisq,yfit=fit)
     resid_p[1:midr-1,i]=mn[1:midr-1,0]*cadence-fit
     fits_p[i,0:1]=res
     fits_p[i,2:3]=sig
     fits_p[i,4]=chisq


     lagstf=reform(lags_ret[*,i,*,0])
     corval=reform(lags_ret[*,i,*,1])
     ;good_p=where(corval gt 0.4) ;Where cross-correlation values less than 0.4  

     mn=(moment(lagstf,dimension=2))[*,0:1] 
     res=poly_fit(dist[1:midr-1],mn[1:midr-1,0]*cadence,1,measure=sqrt(mn[1:midr-1,1])*cadence,sigma=sig,chisq=chisq,yfit=fit)
     resid_r[1:midr-1,i]=mn[1:midr-1,0]*cadence-fit
     fits_r[i,0:1]=res
     fits_r[i,2:3]=sig
     fits_r[i,4]=chisq


ENDFOR

cph_in=1./mean(fits_p[*,1])   ; 1./ mean lag
cph_in_err=(mean(fits_p[*,3])/mean(fits_p[*,1])^2) ; mean error on lag
cph_out=1./mean(fits_r[*,1])
cph_out_err=(mean(fits_r[*,3])/mean(fits_r[*,1])^2)

cph=[[cph_in,cph_out],[cph_in_err,cph_out_err]]

print,'Propagation speed in', cph_in,cph_in_err
print,'Propagation speed out', cph_out,cph_out_err, string(11b)


END



;#########################################################################################
;#########################################################################################
;Main code

PRO loop_analyse_v2_mc,ybast=ybast, tbast=tbast,time_split=time_split,foot=foot,slice=slice,n_slice=n_slice

COMMON loop_param,nx,nz,nt,cadence,pix,mid,sc,ns,mid_eo
COMMON fourier_param,tbas,ybas,half,freq

;inpath='analysis/comp/wave_tracking_output/20120410/'
;outpath=inpath+'test/'
;loop=lp_read(inpath+'SPLINEcuts/velocity/rotTDR_intenscube_V1.fcube')
;restore,inpath+'20120410_loopoutput.sav'
;loop=fltarr(20,320,94)
;for i=0,93 do loop[*,*,i]=reform(output[*,i,1,*])
;fpt=[5,81];fpt=[0,93];


;inpath='analysis/comp/wave_tracking_output/20130914/'
;outpath=inpath
;restore,inpath+'loop_profiles.sav'
;loop=out_v[*,5:170,*]
;loop=rebin(loop[*,0:165,*],17,166/2,95)

inpath='analysis/comp/wave_tracking_output/20130502/'
outpath=inpath
restore,inpath+'loop_profiles.sav'
loop=loop[*,0:175,*]
loop=rebin(loop[*,0:175,*],21,176/2,156/2)


;inpath='analysis/comp/2007/'
;outpath=inpath
;restore,inpath+'loop.sav'




;Phase speed in Mm/s - calculated later in code
;Used for over-plotting phase speeds on w-k plots
cph_in=0.420 ;0.689
cph_out=0.386 ;0.690

;Second half loop
cph_in_2nd=0.391 ; 0.664
cph_out_2nd=0.366 ;0.669


;##########
;Everything below here should be general

set_plot,'ps'
!p.font=0
device,helvetica=1




resolve_routine,'twod_fft',/is_func

;Standard CoMP quantities
cadence=30.     ;s
pix=4.46*0.725  ;Mm
freq_cut=0.005

sz=size(loop)
nx=sz(1) & nt=sz(2) & nz=sz(3)
;Data should be [slice,time,distance]

;Make all series of even length in time - easier!!
IF nt mod 2 EQ 1 THEN loop=loop[0:-1,0:-2,0:-1]
IF nt mod 2 EQ 1 THEN nt=nt-1
;Makes sure split time-series are even too
IF keyword_set(time_split) THEN IF nt/2 MOD 2 EQ 1 THEN BEGIN
loop=loop[0:-1,0:-3,0:-1]
nt=nt-2
ENDIF

;Padding of TD diagrams 
IF n_elements(ybast) eq 0 THEN ybas=0 ELSE ybas=ybast
IF n_elements(tbast) eq 0 THEN tbas=0 ELSE tbas=tbast


IF keyword_set(time_split) THEN nt=nt/2 
IF n_elements(foot) eq 0 THEN fpt=[0,nz-1] ELSE fpt=[foot[0],foot[1]]
IF n_elements(slice) eq 0 THEN sc=[nx/2] ELSE sc=slice
IF n_elements(n_slice) eq 0 THEN n_slice=nx 
ns=[sc-n_slice/2,sc+n_slice/2]


half=(nt-1)/2
X = FINDGEN(half) + 1
IF (nt MOD 2) EQ 0 THEN $
  freq = [0.0, X, nt/2, -nt/2 + X]/(nt*cadence)



;##############################################
;2D FFT
;First half loop

mid=nz ;/2
stop
;Returns average k-w diagrams and variance
kw=loop_pow_mc(loop[ns[0]:ns[1],0:nt-1,fpt[0]:mid-1],1,cph_in,cph_out,0,outpath,inpow,outpow,ratio_mc,freq_n,k_n,halfk)

stop
;Second half of time-series
IF keyword_set(time_split) THEN kw2=loop_pow_mc(loop[ns[0]:ns[1],nt/2:nt*2-1,fpt[0]:mid-1],pow_num_mc,2,cph_in,cph_out,0,outpath,inpow2,outpow2,ratio_mc2)

;##################
;Second half loop

IF nz mod 2 EQ 0 THEN mid_eo=mid ELSE mid_eo=mid+1
kw3=loop_pow_mc(loop[ns[0]:ns[1],0:nt-1,mid_eo:fpt[1]],pow_num_mc,1,cph_in_2nd,cph_out_2nd,1,outpath,inpow3,outpow3,ratio_mc3)

;Second half of time-series
IF keyword_set(time_split) THEN kw4=loop_pow_mc(loop[ns[0]:ns[1],nt/2:nt*2-1,mid_eo:fpt[1]],pow_num_mc,2,cph_in_2nd,cph_out_2nd,1,outpath,inpow4,outpow4,ratio_mc4)



;##############################################
;Plot average power spectra

IF n_elements(kw[*,0,0]) mod 2 EQ 0 THEN midt=[half+2,2*half+1] ELSE midt=[half+1,2*half]
IF n_elements(kw[0,*,0]) mod 2 EQ 0 THEN midk=[halfk+2,2*halfk+1] ELSE midt=[halfk+1,2*halfk]

IF keyword_set(time_split) THEN BEGIN
          aver=(kw[*,*,0]+kw2[*,*,0]+kw3[*,*,0]+kw4[*,*,0])/4.
          sd=sqrt(kw[*,*,1]+kw2[*,*,1]+kw3[*,*,1]+kw4[*,*,1])

ENDIF ELSE BEGIN
          aver=(kw[*,*,0]+kw3[*,*,0])/2.
          sd=sqrt(kw[*,*,1]+kw3[*,*,1])
ENDELSE

averpow=[[[aver]],[[sd]]]


loadct,/silent,13
device,/encapsul,/color,filename=outpath+'power_average.eps'

tvim,alog10(averpow[*,*,0])<(-2.7)>(-4.7),yrange=[0,k_n(halfk)],xrange=[-freq_n[half],freq_n[half]]*1e3,aspect=2.5,$
                          xtitle='Frequency (mHz)',ytitle=' Wavenumber (Mm!e-1!n)'
ks=reverse(freq_n(0:half+1))/cph_in
ks_out=reverse(freq_n(0:half+1))/cph_out  
loadct,0,/silent
oplot,-reverse(freq_n(1:half))*1e3,ks(1:half),color=255,thick=2
oplot,freq_n(1:half)*1e3,reverse(ks_out(1:half)),color=255,thick=2

device,/close



;############################
;Fitting of damping profile




IF keyword_set(time_split) THEN divd=4. ELSE divd=2.
   
IF keyword_set(time_split) THEN BEGIN
   ratio_all=[[[ratio_mc]],[[ratio_mc3]],[[ratio_mc2]],[[ratio_mc4]]] 
ENDIF ELSE ratio_all=[[[ratio_mc]],[[ratio_mc3]]]



nslice=ns[1]-ns[0]+1


;##############################################
;Fitting of damping rate
temp=(size(ratio_all))[1]
sub=temp mod 2
ratio_temp=rebin(ratio_all[0:temp-sub-1,*,*],(temp-sub)/2.,2,divd)

freq_temp=rebin(freq_n(0:temp-sub-1),(temp-sub)/2)

num=where(freq_temp(1:half/2-1) LT freq_cut); and freq_n(1:half) GT 0.001)
mn=min(num)
mx=max(num)
weights=fltarr(mx-mn+1)
weights[*]=1.

num2=where(freq_n(1:half) LT freq_cut); and freq_n(1:half) GT 0.001)
mn2=min(num2)
mx2=max(num2)

reso=[0,0]
reso_pp=[0,0]


FOR i=0,divd-1 DO BEGIN

    IF i mod 2 EQ 0 THEN print,'First leg' ELSE print, 'Second leg'
    IF divd GT 2 THEN IF i LE 2 THEN print,'First half' ELSE print,'Second half'

    res=mpfitfun('myexp',freq_temp(mn+1:mx+1),ratio_temp(mn+1:mx+1,0,i),sqrt(ratio_temp(mn+1:mx+1,1,i)), $
                                                           [1,200.],/quiet,perror=perror,dof=dof,bestnorm=best)
    print,'Exponential fit',res
    print,'Exponential fit error',perror
    print,'Exponential fit chi and DOF',best,dof, string(11b)
    reso=[[reso],[res]]

    res_pp=mpfitfun('pascoe_prof',freq_temp(mn+1:mx+1),ratio_temp(mn+1:mx+1,0,i),sqrt(ratio_temp(mn+1:mx+1,1,i)), $
                                                  [1,200.],/quiet,perror=perror_pp,dof=dof_pp,bestnorm=best_pp)
    print,'Pascoe profile fit',res_pp
    print, 'Pascoe profile fit error', perror_pp
    print, 'Pascoe profile fit chi and DOF', best_pp,dof_pp, string(11b)
    reso_pp=[[reso_pp],[res_pp]]

ENDFOR

    print,'Average fit'
    
    av_rat=mean(ratio_all(mn2+1:mx2+1,0,*),dim=3)
    av_rat_err=sqrt(total(ratio_all(mn2+1:mx2+1,1,*),3))/divd
    res_a=mpfitfun('myexp',freq_n(mn2+1:mx2+1),av_rat,av_rat_err, $
                                                           [1,200.],/quiet,perror=perror,dof=dof,bestnorm=best)
    print,'Exponential fit',res_a
    print,'Exponential fit error',perror
    print,'Exponential fit chi and DOF',best,dof, string(11b)
    

    res_pp_a=mpfitfun('pascoe_prof',freq_n(mn2+1:mx2+1),av_rat,av_rat_err, $
                                                  [1,200.],/quiet,perror=perror_pp,dof=dof_pp,bestnorm=best_pp)
    print,'Pascoe profile fit',res_pp_a
    print, 'Pascoe profile fit error', perror_pp
    print, 'Pascoe profile fit chi and DOF', best_pp,dof_pp, string(11b)
    

verth_mod=res_a[0]*exp(res_a[1]*freq_temp(mn+1:mx+1))
pascoe_mod=res_pp_a[0]*1./(erf(2.*freq_temp(mn+1:mx+1)*res_pp_a[1])/erf(freq_temp(mn+1:mx+1)*res_pp_a[1]) -1.)


loadct,13,/silent
device,/encapsul,/color,filename=outpath+'power_ratio.eps'

!p.multi=0
plot,freq_n(1:half)*1e3,ratio_all[*,0,0],xst=1,xtitle='Frequency (mHz)',ytitle='Power ratio outward to inward',thick=2,/nodata,yrange=[0,max(ratio_all[*,0,*])+1.],charsize=0.8
oploterr,freq_n(1:half)*1e3,ratio_all[*,0,0],sqrt(ratio_all[*,1,0]),thick=2,psym=1
oplot,freq_temp(mn+1:mx+1)*1e3,verth_mod,linestyle=2,color=255,thick=3
oplot,freq_temp(mn+1:mx+1)*1e3,pascoe_mod,linestyle=2,color=150,thick=3
oplot,freq_n(1:half)*1e3,ratio_all[*,0,1],linestyle=2,thick=1
;oploterr,freq_n(mn2+1:mx2+1)*1e3,av_rat,av_rat_err,psym=1,thick=2


IF keyword_set(time_split) THEN BEGIN
;oplot,freq_n(1:half)*1e3,ratio_all[*,0,2],linestyle=3,thick=1
;oplot,freq_n(1:half)*1e3,ratio_all[*,0,3],linestyle=4,thick=1
ENDIF

;plot,freq_n(mn+1:mx+1)*1e3,ratio_all[*,0,0]-verth_mod,linestyle=2,color=255,thick=3,xtitle='Frequency (mHz)',ytitle='Residuals'
;oplot,freq_n(mn+1:mx+1)*1e3,ratio_all[*,0,0]-pascoe_mod,linestyle=2,color=150,thick=3

;Finds probability the two distributions of residuals
;come from same distribution
;i.e., high probability suggests cannot tell models apart
  
device,/close

device,/encapsul,/color,filename=outpath+'ks_test.eps'
kstwo,ratio_all[*,0,0]-verth_mod,ratio_all[*,0,0]-pascoe_mod,d,prob,/plot
print,'Two tailed KS test probability',prob
device,/close


device,/encapsul,/color,filename=outpath+'power_in_vs_out.eps'

!p.multi=[0,2,2]
plot,freq_n(1:half)*1e3,inpow[*,0],xrange=[freq_n(1),freq_n(half)]*1e3,xst=1,/ylog,yrange=[1e-5,1e-2],xtitle='Frequency (mHz)'
oploterr,freq_n(1:half)*1e3,inpow[*,0],sqrt(reverse(inpow[*,1]))
oplot,freq_n(1:half)*1e3,outpow[*,0],linestyle=1,color=255

plot,freq_n(1:half)*1e3,ratio_all[*,0,0],xrange=[freq_n(1),freq_n(half)]*1e3,xst=1,xtitle='Frequency (mHz)'
oploterr,freq_n(1:half)*1e3,ratio_all[*,0,0],sqrt(ratio_all(mn+1:mx+1,1,0))
oplot,freq_n(1:half)*1e3,ratio_all[*,0,1],linestyle=2,color=150


plot,freq_n(1:half)*1e3,(inpow[*,0]),xrange=[freq_n(1),freq_n(half)]*1e3,xst=1,/ylog,yrange=[1e-5,1e-2],xtitle='Frequency (mHz)'
oplot,freq_n(1:half)*1e3,(inpow3[*,0]),color=255
oplot,freq_n(1:half)*1e3,(outpow[*,0]),linestyle=2,color=255
oplot,freq_n(1:half)*1e3,(outpow3[*,0]),linestyle=2

;plot,freq_n(1:half),reverse(inpow[*,0]/outpow3[*,0]),xrange=[freq_n(1),freq_n(half)],xst=1
;oplot,freq_n(1:half),reverse(inpow3[*,0]/outpow[*,0]),color=255


device,/close
!p.multi=0

stop
;##################################
;
;Calculate phase speeds
;Uses Bootstrap method

IF nt mod 2 EQ 0 THEN nt=nt-1 ;deals with even length tracks


num_mc=500 ; number of MC bootstraps to do

print,'Calculating phase speed for first loop leg'
wave_vel_calc,loop,num_mc,[fpt[0],mid-1],n_slice,lags_pro=lags_pro,$
                  lags_ret=lags_ret,fits_p=fits_p,fits_r=fits_r,cph=cph,resid_r=resid_r,resid_p=resid_p


device,/encapsul,/color,filename=outpath+'lag_example.eps'

;choose random plot
cin=fix(randomu(seed)*num_mc)

midr=mid-fpt[0]
dist=(findgen(midr)-midr/2)*pix
IF n_slice mod 2 eq 0 THEN cslic=n_slice+1 ELSE cslic=n_slice

loadct,13,/silent
toplot=reform(lags_pro[*,cin,*,0])


plot,dist,toplot[*,0]*cadence,xtitle='Distance (Mm)',xst=1,ytitle='Lag (s)',yrange=[-120,120]
FOR i=1,cslic-1 DO oplot,dist,toplot[*,i]*cadence
oplot,dist,mean(toplot[*,*],dimension=2)*cadence,thick=3,color=255
oplot,dist,fits_p[cin,0]+fits_p[cin,1]*dist,color=150,linestyle=2
FOR i=0,cslic-1 DO BEGIN
       corval=reform(lags_pro[*,cin,i,1])
       bad_val=where(corval lt 0.4)
       IF n_elements(bad_val) GE 1 THEN oplot,dist[bad_val],toplot[bad_val,i]*cadence,psym=4
ENDFOR
device,/close

stop

print,'Calculating phase speed for second loop leg'
wave_vel_calc,loop,num_mc,[mid_eo,fpt[1]],n_slice,/second_leg,lags_pro=lags_pro2,$
                  lags_ret=lags_ret2,fits_p=fits_p2,fits_r=fits_r2,cph=cph2,resid_r=resid_r2,resid_p=resid_p2


device,/encapsul,/color,filename=outpath+'prop_speed_hist.eps'
!p.multi=[0,2,2]

plothist,1./fits_p[*,1]*1000.,bin=2,xhist,yhist,xtitle='Propagation speed (km s!e-1!n)',charsize=0.85
oplot,xhist,gaussfit(xhist,yhist,nterms=3),linestyle=2,color=255
plots,[cph[0,0],cph[0,0]]*1000.,[0,max(yhist)]

plothist,1./fits_r[*,1]*1000.,bin=2,xhist,yhist,xtitle='Propagation speed (km s!e-1!n)',charsize=0.85
oplot,xhist,gaussfit(xhist,yhist,nterms=3),linestyle=2,color=255
plots,[cph[1,0],cph[1,0]]*1000.,[0,max(yhist)]

loadct,13,/silent
plothist,1./fits_p2[*,1]*1000.,bin=2,xhist,yhist,xtitle='Propagation speed (km s!e-1!n)',charsize=0.85
oplot,xhist,gaussfit(xhist,yhist,nterms=3),linestyle=2,color=255
plots,[cph2[0,0],cph2[0,0]]*1000.,[0,max(yhist)]

plothist,1./fits_r2[*,1]*1000.,bin=2,xhist,yhist,xtitle='Propagation speed (km s!e-1!n)',charsize=0.85
oplot,xhist,gaussfit(xhist,yhist,nterms=3),linestyle=2,color=255
plots,[cph2[1,0],cph2[1,0]]*1000.,[0,max(yhist)]



device,/close
!p.multi=0
set_plot,'x'



END










