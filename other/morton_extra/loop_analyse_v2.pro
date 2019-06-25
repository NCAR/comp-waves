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
FUNCTION loop_pow,data,times,pos,tbas,ybas,loop_num,cph_in,cph_out,lgno,outpath,inpow,outpow,$
                  freq_n,k_n,half,halfk

COMMON loop_param,nx,nz,nt,cadence,pix,mid,sc,ns,mid_eo

IF lgno EQ 1 THEN rot=4 ELSE rot=1 


;set up array
in=transpose( reform(data[sc,times[0]:times[1],pos[0]:pos[1]]) ) 
im=TWOD_fft_comp(in,cadence,pix,freq_n=freq_n,k_n=k_n,tbas=tbas,ybas=ybas,/han)



sz=size(im)
nyt=sz[1] & ntt=sz[2]
half=(ntt-1)/2
halfk=(nyt-1)/2

IF ntt mod 2 eq 0 THEN lbound=half+2 ELSE lbound=half+1
IF ntt mod 2 eq 0 THEN sub=1 ELSE sub=0

nslice=ns[1]-ns[0]+1
twod_pow=fltarr((nyt-1)/2+1,ntt-sub,nslice) ;create odd length array if ntt even

midr=(size(in))[1]
FOR i=ns[0],ns[1] DO BEGIN

    in=transpose(reform(data[i,times[0]:times[1],pos[0]:pos[1]])) 
    FOR j=0,midr-1 DO in[j,*]=in[j,*]-mean(in[j,*])
    FOR j=0,nt/2-1 DO in[*,j]=in[*,j]-mean(in[*,j])
    
    im=TWOD_fft_comp(in,cadence,pix,tbas=tbas,ybas=ybas,/han)
    
    twod_pow[1:halfk,0:half-1,i-nslice/2]=2.*reverse(im[1:halfk,1:half],2)
    twod_pow[1:halfk,half+1:ntt-1-sub,i-nslice/2]=2.*reverse(im[1:halfk,lbound:ntt-1],2)
     
   
ENDFOR

;Average k-w power and variance from all slices
kw=(moment(twod_pow,dim=3))[*,*,0:1]
kw=[[[rotate(kw[*,*,0],rot)]],[[rotate(kw[*,*,1],rot)]]]

;Interpolate power along w/k=c_ph

ks=reverse(freq_n(0:half+1))/cph_in
ks_out=reverse(freq_n(0:half+1))/cph_out  

inpow=interpolate(kw[*,*,0],freq_n(0:half-1)*cadence*ntt,ks(1:half)*pix*nyt)
inpow_var=interpolate(kw[*,*,1],freq_n(0:half-1)*cadence*ntt,ks(1:half)*pix*nyt)
outpow=interpolate(kw[*,*,0],ntt-sub-freq_n(1:half)*cadence*ntt,ks_out(1:half)*pix*nyt)
outpow_var=interpolate(kw[*,*,1],ntt-sub-freq_n(1:half)*cadence*ntt,ks_out(1:half)*pix*nyt)


inpow=[[reform(inpow)],[reform(inpow_var)]]
outpow=[[outpow],[outpow_var]]


;##############################################
;Plot loop power
loadct,13,/silent
device,/encapsul,/color,filename=outpath+'power_'+strcompress(loop_num)+'st_half_leg'+strcompress(lgno)+'.eps'

vals=moment(alog10(kw[*,*,0])>(-10)<0.)
bot=-6. ;vals[0]-sqrt(vals[1])
top=-2.5 ;vals[0]+sqrt(vals[1])

tvim,alog10(kw[*,*,0])<(top)>(bot),yrange=[0,k_n(halfk)],xrange=[-freq_n[half],freq_n[half]],aspect=2.5,$
                          xtitle='Frequency (Hz)',ytitle=' Wavenumber (Mm!e-1!n)'

loadct,0,/silent
oplot,-reverse(freq_n(1:half)),ks(1:half),color=255,thick=3
oplot,freq_n(1:half),reverse(ks_out(1:half)),color=255,thick=3

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
n = midr  ;number of points along track for initial guess (make odd)



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

PRO loop_analyse_v2,ybas=ybas, tbas=tbas,time_split=time_split,foot=foot,slice=slice,n_slice=n_slice

inpath='data/comp/wave_tracking_output/20120410/'
outpath=inpath+'test/'
;loop=lp_read(inpath+'SPLINEcuts/velocity/rotTDR_intenscube_V1.fcube')
restore,inpath+'20120410_loopoutput.sav'
loop=fltarr(20,320,94)
for i=0,93 do loop[*,*,i]=reform(output[*,i,1,*])
fpt=[0,93];fpt=[5,81]


;inpath='data/comp/wave_tracking_output/20130914/'
;outpath=inpath
;restore,inpath+'loop_profiles.sav'
;loop=out_v[*,0:170,*]



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


COMMON loop_param,nx,nz,nt,cadence,pix,mid,sc,ns,mid_eo

;Standard CoMP quantities
cadence=30.     ;s
pix=4.46*0.725 ;Mm

sz=size(loop)
nx=sz(1) & nt=sz(2) & nz=sz(3)



IF n_elements(ybas) eq 0 THEN ybas=0
IF n_elements(tbas) eq 0 THEN tbas=0

IF keyword_set(time_split) THEN ts=[nt/2-1,nt-1] ELSE ts=nt-1
IF n_elements(foot) eq 0 THEN fpt=[0,nz-1] ELSE fpt=[foot[0],foot[1]]
IF n_elements(slice) eq 0 THEN sc=[nx/2] ELSE sc=slice
IF n_elements(n_slice) eq 0 THEN n_slice=10 
ns=[sc-n_slice/2,sc+n_slice/2]

stop
; define temporal apodization
apodt = fltarr(nt)+1
apod=0.1
apodrimt = nt*apod
apodt[0] = (sin(!pi/2.*findgen(apodrimt)/apodrimt))^2
apodt = apodt*shift(rotate(apodt,2),1) 
CPGT=total(apodt)/n_elements(apodt) ;Coherent Power Gain - correction factor needed for power after apodisation



half=(nt-1)/2
X = FINDGEN(half) + 1
IF (nt MOD 2) EQ 0 THEN $
  freq = [0.0, X, nt/2, -nt/2 + X]/(nt*cadence) $
ELSE $
  freq = [0.0, X, -(nt/2 + 1) + X]/(nt*cadence)



;##############################################
;2D FFT
;First half loop

mid=nz/2

print,tbas

;Returns average k-w diagrams and variance
kw=loop_pow(loop,[0,ts[0]],[fpt[0],mid-1],tbas,ybas,1,cph_in,cph_out,0,outpath,inpow,outpow,freq_n,k_n,half,halfk)

;Second half of time-series
IF keyword_set(time_split) THEN kw2=loop_pow(loop,[ts[0]+1,ts[1]],[fpt[0],mid-1],tbas,ybas,2,cph_in,cph_out,0,outpath,inpow2,outpow2)


;##################
;Second half loop

IF nz mod 2 EQ 0 THEN mid_eo=mid ELSE mid_eo=mid+1
kw3=loop_pow(loop,[0,ts[0]],[mid_eo,fpt[1]],tbas,ybas,1,cph_in_2nd,cph_out_2nd,1,outpath,inpow3,outpow3)

;Second half of time-series
IF keyword_set(time_split) THEN kw4=loop_pow(loop,[ts[0]+1,ts[1]],[mid,fpt[1]],tbas,ybas,2,cph_in_2nd,cph_out_2nd,1,outpath,inpow4,outpow4)



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

tvim,alog10(averpow[*,*,0])<(-2.5)>(-4.5),yrange=[0,k_n(halfk)],xrange=[-freq_n[half],freq_n[half]],aspect=2.5,$
                          xtitle='Frequency (Hz)',ytitle=' Wavenumber (Mm!e-1!n)'
ks=reverse(freq_n(0:half+1))/cph_in
ks_out=reverse(freq_n(0:half+1))/cph_out  
loadct,0,/silent
oplot,-reverse(freq_n(1:half)),ks(1:half)
oplot,freq_n(1:half),reverse(ks_out(1:half))

device,/close



;############################
;Fitting of damping profile



;All power below 0.05 Mm-1, i.e. as in Verth et al.,
in=where(k_n lt 0.05)
IF keyword_set(time_split) THEN BEGIN
  powerin=[[total(kw[0:half,1:max(in),0],2)],[total(kw2[0:half,1:max(in),0],2)],[total(kw3[0:half,1:max(in),0],2)],[total(kw4[0:half,1:max(in),0],2)]] 
  powerin_var=[[total(kw[0:half,1:max(in),1],2)],[total(kw2[0:half,1:max(in),1],2)],[total(kw3[0:half,1:max(in),1],2)],[total(kw4[0:half,1:max(in),1],2)]] 

  powerout=[[total(kw[midt[0]:midt[1],1:max(in),0],2)],[total(kw2[midt[0]:midt[1],1:max(in),0],2)],[total(kw3[midt[0]:midt[1],1:max(in),0],2)],[total(kw4[midt[0]:midt[1],1:max(in),0],2)]]
 powerout_var=[[total(kw[midt[0]:midt[1],1:max(in),1],2)],[total(kw2[midt[0]:midt[1],1:max(in),1],2)],[total(kw3[midt[0]:midt[1],1:max(in),1],2)],[total(kw4[midt[0]:midt[1],1:max(in),1],2)]]

ENDIF ELSE BEGIN

    powerin=[[total(kw[0:half,1:max(in),0],2)],[total(kw3[0:half,1:max(in),0],2)]]
    powerin_var=[[total(kw[0:half,1:max(in),1],2)],[total(kw3[0:half,1:max(in),1],2)]] 
    powerout=[[total(kw[midt[0]:midt[1],1:max(in),0],2)],[total(kw3[midt[0]:midt[1],1:max(in),0],2)]]
    powerout_var=[[total(kw[midt[0]:midt[1],1:max(in),1],2)],[total(kw3[midt[0]:midt[1],1:max(in),1],2)]]
ENDELSE

IF keyword_set(time_split) THEN divd=4. ELSE divd=2.

powerin=mean(powerin,dim=2)
powerin_var=total(powerin_var,2)/divd^2

powerout=mean(powerout,dim=2)
powerout_var=total(powerout_var,2)/divd^2

ratio_all=reverse(powerin)/powerout
err_all=ratio_all*reverse(sqrt(reverse(powerin_var/powerin^2)+powerout_var/powerout^2))



;Power along cph lines
IF keyword_set(time_split) THEN BEGIN
    totin= [[inpow[*,0]],[inpow2[*,0]],[inpow3[*,0]],[inpow4[*,0]]]
    totin_var= [[inpow[*,1]],[inpow2[*,1]],[inpow3[*,1]],[inpow4[*,1]]]
    totout= [[outpow[*,0]],[outpow2[*,0]],[outpow3[*,0]],[outpow4[*,0]]]
    totout_var= [[outpow[*,1]],[outpow2[*,1]],[outpow3[*,1]],[outpow4[*,1]]]

ENDIF ELSE BEGIN
    totin=[[inpow[*,0]],[inpow3[*,0]]] ;mean values
    totin_var=[[inpow[*,1]],[inpow3[*,1]]] ; variance values
    totout= [[outpow[*,0]],[outpow3[*,0]]]
    totout_var= [[outpow[*,1]],[outpow3[*,1]]]
ENDELSE





; mean & SD of all slices and foot-points
totin=[[mean(totin,dim=2)],[sqrt(total(totin_var,2))/divd]] 
totout=[[mean(totout,dim=2)],[sqrt(total(totout_var,2))/divd]]



ratio=reverse(totin[*,0]/totout[*,0])
err=ratio*reverse(sqrt(totin[*,1]^2/totin[*,0]^2+totout[*,1]^2/totout[*,0]^2))

num=where(freq_n(1:half) lt 0.007)
weight=fltarr(max(num)+1)
weight[*]=[1.]

;Not sure we can trust calculated variances
res=mpfitfun('myexp',freq_n(1:max(num)+1),ratio(0:max(num)),err(0:max(num)),[1,200.],/quiet,perror=perror,dof=dof,bestnorm=best,weight=weight)
print,'Exponential fit',res
print,'Exponential fit error',perror
print,'Exponential fit chi and DOF',best,dof, string(11b)


res_pp=mpfitfun('pascoe_prof',freq_n(1:max(num)+1),ratio(0:max(num)),err(0:max(num)),[1,200.],/quiet,perror=perror_pp,dof=dof_pp,bestnorm=best_pp,weight=weight)
print,'Pascoe profile fit',res_pp
print, 'Pascoe profile fit error', perror_pp
print, 'Pascoe profile fit chi and DOF', best_pp,dof_pp, string(11b)

loadct,13,/silent
device,/encapsul,/color,filename=outpath+'power_ratio.eps'

plot,freq_n(1:half),ratio,xst=1,xtitle='Frequency (Hz)',ytitle='Power ratio outward to inward',thick=3
;oploterr,freq_n(1:half),ratio,err,thick=3
oplot,freq_n(1:max(num)),res[0]*exp(res[1]*freq_n(1:max(num))),linestyle=2,color=255,thick=3
;oploterr,freq_n(1:half),ratio_all,err_all,color=150,thick=3,errcolor=150
oplot,freq_n(1:max(num)),res_pp[0]*1./(erf(2.*freq_n(1:max(num))*res_pp[1])/erf(freq_n(1:max(num))*res_pp[1]) -1.),linestyle=2,color=150,thick=3


device,/close

device,/encapsul,/color,filename=outpath+'power_in_vs_out.eps'

!p.multi=[0,2,2]
plot,freq_n(1:half),reverse(powerin[*]),xrange=[freq_n(1),freq_n(half)],xst=1,/ylog,yrange=[1e-5,1]
oploterr,freq_n(1:half),reverse(powerin[*]),sqrt(reverse(powerin_var[*]))
oplot,freq_n(1:half),powerout[*],linestyle=1,color=255

plot,freq_n(1:half),ratio_all,xrange=[freq_n(1),freq_n(half)],xst=1
oploterr,freq_n(1:half),ratio_all,err_all
oplot,freq_n(1:half),ratio,color=150
device,/close



!p.multi=0


;##################################
;
;Calculate phase speeds
;Uses Bootstrap method

IF nt mod 2 EQ 0 THEN nt=nt-1 ;deals with even length tracks


num_mc=200 ; number of MC bootstraps to do

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










;++++++++++++++++++++++++++++++++++++++++++++++++++

;OLD CODE replaced by wave_vel_calc
;randerr=randomn(systime,nt,mid,num_mc) ;generate random errors
;lags_pro=fltarr(mid,num_mc,n_slice+1)
;lags_ret=fltarr(mid,num_mc,n_slice+1)
;amplitude2=fltarr(mid,n_slice+1,2)
;NEEDED FOR RECALCULATING AMPLITUDE
;half=(nt-1)/2
;X = FINDGEN(half) + 1
;IF (nt MOD 2) EQ 0 THEN $
;  f = [0.0, X, nt/2, -nt/2 + X]/(nt*cadence) $
;ELSE $
;  f = [0.0, X, -(nt/2 + 1) + X]/(nt*cadence)

;fin=where(f gt 0 and f lt 0.007)
;vph=0.390
;xi=4.2759
;fac=pix/vph/xi


;FOR k=ns[0],ns[1] DO BEGIN

;    vel=rotate(reform(loop[k,0:nt-1,mid_eo:fpt[1]]),7)
;    FOR i=0,mid-1 DO vel[*,i]=vel[*,i]-mean(vel[*,i])

;    trans=fft(vel)
;    pro_trans=trans             ;select prograde waves (assume nt even, npt odd)
;    pro_trans(1:nt/2-1,0:mid/2)=0.
;    pro_trans(nt/2+1:nt-1,mid/2+1:mid-1)=0.
;    pro_vel=float(fft(pro_trans,/inverse))
;    err=fltarr(mid)
;    FOR i=0,mid-1 DO err[i]=sqrt(mean((pro_vel[*,i]-smooth(pro_vel[*,i],3,/edge_truncate))^2))
;    amplitude2[*,k-ns[0]]=sqrt(mean(smooth(pro_vel,[3,1])^2,dim=1))

;    errs=randerr*transpose(rebin(err,mid,nt,num_mc),[1,0,2])
;    lags_pro[*,*,k-ns[0]]=calc_speed(n,mid,nt,icent,pro_vel,errs,num_mc)
    
 
;    ret_trans=trans             ;select retrograde waves
;    ret_trans(1:nt/2-1,mid/2+1:mid-1)=0.
;    ret_trans(nt/2+1:nt-1,0:mid/2)=0.
;    ret_vel=float(fft(ret_trans,/inverse))
;    FOR i=0,mid-1 DO err[i]=sqrt(mean((ret_vel[*,i]-smooth(ret_vel[*,i],3,/edge_truncate))^2))
    
;    errs=randerr*transpose(rebin(err,mid,nt,num_mc),[1,0,2])
;    lags_ret[*,*,k-ns[0]]=calc_speed(n,mid,nt,icent,ret_vel,errs,num_mc) 

    ;For recalculating amplitude using damping measurements
;    amplitude2[*,k-ns[0],0]=sqrt(mean(smooth(pro_vel,[3,1])^2,dim=1))
;    FOR j=0,mid-1 DO BEGIN
;      dum=abs(fft(smooth(pro_vel[*,j],3)))^2
;      amplitude2[j,k-ns[0],1]=sqrt( total( 2.*dum[fin]*exp(2.*fac*j*f[fin]) ) )
;    ENDFOR


;ENDFOR
