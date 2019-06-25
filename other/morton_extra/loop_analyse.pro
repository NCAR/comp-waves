;LOOP ANALYSE.pro
;
;PURPOSE: Used to analyse coronal loops from CoMP data. Will separate inward and outward
;propagating power for both loop legs to enable estimate of damping & propagation speeds.
;
;
;
;
;

;Calculates the propagation of velocity signals via cross-correlation
;Uses noise estimated from signal to produce bootstrap repetitions of
;CC for each row.
;Borrowed and modified from compute_speed.pro by S. Tomczyk
FUNCTION calc_speed, n,mid,nt,icent,pro_vel,errs,num_mc

;##############################################################
;Cross-correlate adjacent rows in time-distance map
;##############################################################
cent=fltarr(n)
FOR i=-n/2,n/2 DO BEGIN
  ccor=ccpeak(pro_vel(*,icent),pro_vel(*,icent+i))
  cent(i+n/2)=ccor
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
dcent=fltarr(mid)+.1

FOR j=0,num_mc-1 DO BEGIN
    FOR i=0,mid-1 DO BEGIN
        ccor=ccpeak(col,smooth(pro_vel(*,i),3,/edge_truncate)+errs(*,i,j))
        cent(i,j)=ccor
  
   ENDFOR
ENDFOR

return,cent
END


;Main code
PRO loop_analyse,ybas=ybas, tbas=tabs,time_split=time_split

inpath='data/comp/wave_tracking_output/20120410/'
outpath=inpath+'damping/'
;loop=lp_read(inpath+'SPLINEcuts/velocity/rotTDR_intenscube_V1.fcube')
restore,inpath+'20120410_loopoutput.sav'

loop=fltarr(20,320,94)

for i=0,93 do loop[*,*,i]=reform(output[*,i,1,*])

fpt=[0,93];fpt=[5,81]

;Phase speed in Mm/s - calculated later in code
;Used for over-plotting phase speeds on w-k plots
cph_in=0.420
cph_out=0.386


;##########
;Everything below here should be general

;Standard CoMP quantities
cadence=30.     ;s
pix=4.46*0.725 ;Mm

sz=size(loop)
nx=sz(1) & nt=sz(2) & nz=sz(3)

IF n_elements(ybas) eq 0 THEN ybas=1
IF n_elements(tbas) eq 0 THEN tbas=1
IF n_elements(time_split) EQ 0 THEN t1=nt-1




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


ny=nz/2-fpt[0]+1
mid=nz/2

;set up array
in=transpose(reform(loop[10,0:159,fpt[0]:mid-1])) 
im=TWOD_fft_comp(in,cadence,pix,freq_n=freq_n,k_n=k_n,tbas=tbas,ybas=ybas,/han)
sz=size(im)
nyt=sz[1] & ntt=sz[2]
half=(ntt-1)/2
halfk=(nyt-1)/2
twod_pow=fltarr((nyt-1)/2+1,ntt,11)

loop_fft=fltarr(ny,nt,nx)

FOR i=5,15 DO BEGIN

    in=transpose(reform(loop[i,0:159,fpt[0]:mid-1])) 
    FOR j=0,mid-1 DO in[j,*]=in[j,*]-mean(in[j,*])
    FOR j=0,nt/2-1 DO in[*,j]=in[*,j]-mean(in[*,j])
    
    im=TWOD_fft_comp(in,cadence,pix,tbas=tbas,ybas=ybas,/han)
    
    twod_pow[1:(nyt-1)/2,0:(ntt-1)/2-1,i-5]=2.*reverse(im[1:(nyt-1)/2,1:(ntt-1)/2],2)
    twod_pow[1:(nyt-1)/2,(ntt-1)/2+1:ntt-1,i-5]=2.*reverse(im[1:(nyt-1)/2,(ntt-1)/2+1:ntt-1],2)
     

    in=transpose(reform(loop[i,*,fpt[0]:mid-1])) 
    FOR j=0,mid-1 DO in[j,*]=in[j,*]-mean(in[j,*])
    FOR j=0,nt/2-1 DO in[*,j]=in[*,j]-mean(in[*,j])
    loop_fft[1:(ny-1)/2,0:(nt-1)/2-1,i-5]=2.*reverse((abs(fft(in))^2)[1:(ny-1)/2,1:(nt-1)/2],2)
    loop_fft[1:(ny-1)/2,(nt-1)/2+1:nt-1,i-5]=2.*reverse((abs(fft(in))^2)[1:(ny-1)/2,(nt-1)/2+1:nt-1],2)
   
ENDFOR
firsthalf=twod_pow
kw=(moment(twod_pow,dim=3))[*,*,0:1]
kw=[[[rotate(kw[*,*,0],1)]],[[rotate(kw[*,*,1],1)]]]
loadct,13,/silent

set_plot,'ps'
!p.font=0
device,helvetica=1

;##############################################
;Plot 1st half loop power 1st time interval
device,/encapsul,/color,filename=outpath+'power_1st_half.eps'

tvim,alog10(kw[*,*,0])<(-2.5)>(-4.5),yrange=[0,k_n(halfk)],xrange=[-freq_n[half+1],freq_n[half]],aspect=2.5,$
                          xtitle='Frequency (Hz)',ytitle=' Wavenumber (Mm!e-1!n)'
ks=reverse(freq_n(0:half+1))/cph_in
ks_out=reverse(freq_n(0:half+1))/cph_out  
loadct,0,/silent
oplot,-reverse(freq_n(0:half+1)),ks,color=255,thick=3
oplot,freq_n(0:half),reverse(ks_out),color=255,thick=3

device,/close

inpow=interpolate(kw[*,*,0],freq_n(1:half)*cadence*ntt,ks(1:half)*pix*nyt)
outpow=interpolate(kw[*,*,0],ntt-freq_n(1:half)*cadence*ntt,ks_out(1:half)*pix*nyt)

fcub=reverse(rebin(freq_n[0:half+1],half+2,halfk+1))
kcub=transpose(rebin(k_n[0:halfk],halfk+1,half+2))
phcub=fcub/kcub
phcub=[phcub,reverse(phcub[1:half,*])]


kwall=fltarr(ntt,halfk+1,11)
FOR i=0,10 DO kwall[*,*,i]=rotate(twod_pow[*,*,i],1)

in1=where(phcub lt 0.3)
in2=where(phcub gt 0.9)
FOR i=0,10 DO BEGIN
    dum=kwall[*,*,i]
    dum[in1]=0
    dum[in2]=0
    kwall[*,*,i]=dum

ENDFOR
pow=total(kwall,2)

twod_pow_ptt=fltarr((nyt-1)/2+1,ntt,11)


;##############################################
;Second half of time-series

FOR i=5,15 DO BEGIN

    in=transpose(reform(loop[i,160:319,fpt[0]:mid-1])) 
    FOR j=0,mid-1 DO in[j,*]=in[j,*]-mean(in[j,*])
    FOR j=0,nt/2-1 DO in[*,j]=in[*,j]-mean(in[*,j])
    ;atrous,in,decomp=a,n_scales=1
    ;in=in-a[*,*,1]
    im=TWOD_fft_comp(in,cadence,pix,tbas=tbas,ybas=ybas,/han)
    
    twod_pow_ptt[1:(nyt-1)/2,0:(ntt-1)/2-1,i-5]=2.*reverse(im[1:(nyt-1)/2,1:(ntt-1)/2],2)
    twod_pow_ptt[1:(nyt-1)/2,ntt/2+1:ntt-1,i-5]=2.*reverse(im[1:(nyt-1)/2,ntt/2+1:ntt-1],2)
    ;stop  
   
ENDFOR
firsthalf2=twod_pow_ptt
kw2=(moment(twod_pow_ptt,dim=3))[*,*,0:1]
kw2=[[[rotate(kw2[*,*,0],1)]],[[rotate(kw2[*,*,1],1)]]]
inpow2=interpolate(kw2[*,*,0],freq_n(1:half)*cadence*ntt,ks(1:half)*pix*nyt)
outpow2=interpolate(kw2[*,*,0],ntt-freq_n(1:half)*cadence*ntt,ks_out(1:half)*pix*nyt)


;##############################################
;2D FFT
;Second half loop
;Phase speed in Mm/s
cph_in=0.391
cph_out=0.366

ny=fpt[1]-(nz/2+1)+1

;set up array
in=transpose(reform(loop[10,0:159,mid:fpt[1]])) 
im=TWOD_fft_comp(in,cadence,pix,freq_n=freq_n,k_n=k_n,tbas=tbas,ybas=ybas,/han)
sz=size(im)
nyt=sz[1] & ntt=sz[2]
half=(ntt-1)/2
halfk=(nyt-1)/2
twod_pow=fltarr((nyt-1)/2+1,ntt,11)

FOR i=5,15 DO BEGIN

    in=transpose(reform(loop[i,0:159,mid:fpt[1]])) 
    FOR j=0,mid-1 DO in[j,*]=in[j,*]-mean(in[j,*])
    FOR j=0,nt/2-1 DO in[*,j]=in[*,j]-mean(in[*,j])
    im=TWOD_fft_comp(in,cadence,pix,tbas=tbas,ybas=ybas,/han)
 
    
    twod_pow[1:(nyt-1)/2,0:(ntt-1)/2-1,i-5]=2.*reverse(im[1:(nyt-1)/2,1:(ntt-1)/2],2)
    twod_pow[1:(nyt-1)/2,(ntt-1)/2+1:ntt-1,i-5]=2.*reverse(im[1:(nyt-1)/2,(ntt-1)/2+1:ntt-1],2)
   
ENDFOR
secondhalf=twod_pow
kw3=(moment(twod_pow,dim=3))[*,*,0:1]
kw3=[[[rotate(kw3[*,*,0],4)]],[[rotate(kw3[*,*,1],4)]]]
loadct,13,/silent



;##############################################
;Plot 2nd half loop power 2nd time interval
device,/encapsul,/color,filename=outpath+'power_2nd_half.eps'

tvim,alog10(kw3[*,*,0])<(-2.5)>(-4.5),yrange=[0,k_n(halfk)],xrange=[-freq_n[half+1],freq_n[half]],aspect=2.5,$
                          xtitle='Frequency (Hz)',ytitle=' Wavenumber (Mm!e-1!n)'
ks=reverse(freq_n(0:half+1))/cph_in
ks_out=reverse(freq_n(0:half+1))/cph_out  
loadct,0,/silent
oplot,-reverse(freq_n(0:half+1)),ks,color=255,thick=3
oplot,freq_n(0:half),reverse(ks_out),color=255,thick=3

device,/close

inpow3=interpolate(kw3[*,*,0],freq_n(1:half)*cadence*ntt,ks(1:half)*pix*nyt)
outpow3=interpolate(kw3[*,*,0],ntt-freq_n(1:half)*cadence*ntt,ks_out(1:half)*pix*nyt)


;Second part of time-series
twod_pow=fltarr((nyt-1)/2+1,ntt,11)

FOR i=5,15 DO BEGIN

    in=transpose(reform(loop[i,160:319,mid:fpt[1]])) 
    FOR j=0,mid-1 DO in[j,*]=in[j,*]-mean(in[j,*])
    FOR j=0,nt/2-1 DO in[*,j]=in[*,j]-mean(in[*,j])
    im=TWOD_fft_comp(in,cadence,pix,tbas=tbas,ybas=ybas,/han)
    
    twod_pow[1:(nyt-1)/2,0:(ntt-1)/2-1,i-5]=2.*reverse(im[1:(nyt-1)/2,1:(ntt-1)/2],2)
    twod_pow[1:(nyt-1)/2,ntt/2+1:ntt-1,i-5]=2.*reverse(im[1:(nyt-1)/2,ntt/2+1:ntt-1],2)
    
   
ENDFOR
secondhalf2=twod_pow
kw4=(moment(twod_pow,dim=3))[*,*,0:1]
kw4=[[[rotate(kw4[*,*,0],4)]],[[rotate(kw4[*,*,1],4)]]]
inpow4=interpolate(kw4[*,*,0],freq_n(1:half)*cadence*ntt,ks(1:half)*pix*nyt)
outpow4=interpolate(kw4[*,*,0],ntt-freq_n(1:half)*cadence*ntt,ks(1:half)*pix*nyt)


loadct,13,/silent

;##############################################
;Plot average power spectra

tot_test=fltarr(256,32,44)
FOR i=0,10 DO tot_test[*,*,i]=rotate(firsthalf[*,*,i],1)
FOR i=0,10 DO tot_test[*,*,i+11]=rotate(firsthalf2[*,*,i],1)
FOR i=0,10 DO tot_test[*,*,i+22]=rotate(secondhalf[*,*,i],4)
FOR i=0,10 DO tot_test[*,*,i+33]=rotate(secondhalf2[*,*,i],4)

device,/encapsul,/color,filename=outpath+'power_average.eps'
averpow=(moment(tot_test,dim=3))[*,*,0:1]
tvim,alog10(averpow[*,*,0])<(-2.5)>(-4.5),yrange=[0,k_n(halfk)],xrange=[-freq_n[half+1],freq_n[half]],aspect=2.5,$
                          xtitle='Frequency (Hz)',ytitle=' Wavenumber (Mm!e-1!n)'
ks=reverse(freq_n(0:half+1))/cph_in
ks_out=reverse(freq_n(0:half+1))/cph_out  
loadct,0,/silent
oplot,-reverse(freq_n(0:half+1)),ks
oplot,freq_n(0:half),reverse(ks_out)

device,/close

stop

avinpow=interpolate(averpow[*,*,0],freq_n(1:half)*cadence*ntt,ks(1:half)*pix*nyt)
avoutpow=interpolate(averpow[*,*,0],ntt-freq_n(1:half)*cadence*ntt,ks(1:half)*pix*nyt)

avinpow_er=interpolate(averpow[*,*,1],freq_n(1:half)*cadence*ntt,ks(1:half)*pix*nyt)
avoutpow_er=interpolate(averpow[*,*,1],ntt-freq_n(1:half)*cadence*ntt,ks(1:half)*pix*nyt)

;############################
;Fitting of damping profile

device,/encapsul,/color,filename=outpath+'power_in_vs_out.eps'
in=where(k_n lt 0.05)
powerin=[[total(kw[0:half,1:max(in),0],2)],[total(kw2[0:half,1:max(in),0],2)],[total(kw3[0:half,1:max(in),0],2)],[total(kw4[0:half,1:max(in),0],2)]] 
powerout=[[total(kw[half+2:2*half+1,1:max(in),0],2)],[total(kw2[half+2:2*half+1,1:max(in),0],2)],[total(kw3[half+2:2*half+1,1:max(in),0],2)],[total(kw4[half+2:2*half+1,1:max(in),0],2)]]
powerin=(moment(powerin,dim=2))[*,0:1]
powerout=(moment(powerout,dim=2))[*,0:1]
;errpow=sqrt( )
loadct,13
plot,freq_n(1:half),reverse(powerin[*,0]),xrange=[freq_n(1),freq_n(half)],xst=1,/ylog,yrange=[1e-5,1]
oploterr,freq_n(1:half),reverse(powerin[*,0]),sqrt(reverse(powerin[*,1]))
oplot,freq_n(1:half),(powerout[*,0]),linestyle=1
oploterr,freq_n(1:half),reverse(avinpow[*]),sqrt(reverse(avinpow_er[*]))
oploterr,freq_n(1:half),reverse(avoutpow[*]),sqrt(reverse(avoutpow_er[*])),color=255
device,/close



totin=[[inpow],[inpow2],[inpow3],[inpow4]] ;;
totout= [[outpow],[outpow2],[outpow3],[outpow4]]
totin=((moment(totin,dim=2))[*,0:1])
totout=(moment(totout,dim=2))[*,0:1]

ratio=reverse(totin[*,0]/totout[*,0])
err=ratio*reverse(sqrt(totin[*,1]/totin[*,0]^2+totout[*,1]/totout[*,0]^2))

num=where(freq_n(1:half) lt 0.006)

res=mpfitfun('myexp',freq_n(1:max(num)+1),ratio(0:max(num)),err(0:max(num)),[1,200.],/quiet,perror=perror,dof=dof,bestnorm=best)
print,res
print,perror
print,best,dof

res_pp=mpfitfun('pascoe_prof',freq_n(1:max(num)+1),ratio(0:max(num)),err(0:max(num)),[1,200.],/quiet,perror=perror_pp,dof=dof_pp,bestnorm=best_pp)
print,res_pp
print,perror_pp
print,best_pp,dof_pp

device,/encapsul,/color,filename=outpath+'power_ratio.eps'

plot,freq_n(1:half),ratio,psym=1,xst=1,xtitle='Frequency (Hz)',ytitle='Power ratio outward to inward'
oploterr,freq_n(1:half),ratio,err,thick=3
oplot,freq_n(1:max(num)),res[0]*exp(res[1]*freq_n(1:max(num))),linestyle=2,color=255,thick=3
;oplot,freq_n(1:half),reverse(avinpow/avoutpow)


oplot,freq_n(1:max(num)),res_pp[0]*1./(erf(2.*freq_n(1:max(num))*res_pp[1])/erf(freq_n(1:max(num))*res_pp[1]) -1.),linestyle=2,color=150,thick=3


device,/close

stop



;The MC Bootstrap version

;ny=nz/2-fpt[0]+1

;set up array
;in=transpose(reform(loop[10,*,fpt[0]:nz/2])) 
;im=TWOD_fft_comp(in,cadence,pix,freq_n=freq_n,k_n=k_n,tbas=2,ybas=3,/han)
;sz=size(im)
;nyt=sz[1] & ntt=sz[2]
;half=(ntt-1)/2
;halfk=(nyt-1)/2
;twod_pow=fltarr((nyt-1)/2+1,ntt,11,500)

;rndser=randomn(systime,ny,nt,11,500)

;FOR i=5,15 DO BEGIN
;    in=transpose(reform(loop[i,*,fpt[0]:nz/2])) 
;    FOR j=0,ny-1 DO in[j,*]=in[j,*]-mean(in[j,*])
;    FOR j=0,nt-1 DO in[*,j]=in[*,j]-mean(in[*,j])
    
;   sds=fltarr(ny)
;   FOR j=0,ny-1 DO sds[j]=sqrt(mean((reform(in[j,*])-smooth(reform(in[j,*]),3,/edge_truncate))^2))   
   
;  FOR k=0,499 DO BEGIN
;    ser=smooth(in,3,/edge_truncate)+reform(rndser[*,*,i-5,k])*rebin(sds,ny,nt)
;    im=TWOD_fft_comp(ser,cadence,pix,tbas=2,ybas=3,/han)
;    twod_pow[1:(nyt-1)/2,0:(ntt-1)/2-1,i-5,k]=2.*reverse(im[1:(nyt-1)/2,1:(ntt-1)/2],2)
;    twod_pow[1:(nyt-1)/2,ntt/2+1:ntt-1,i-5,k]=2.*reverse(im[1:(nyt-1)/2,ntt/2+1:ntt-1],2)
    
;  ENDFOR    

   
;ENDFOR

;kw=rotate(total(twod_pow,3),1)
;loadct,13,/silent
;tvim,alog10(kw)<(-2)>(-4),yrange=[0,k_n((nyt-1)/2-1)],xrange=[-f[0:nt/2],f[0:nt/2]]
;ks=reverse(freq_n)/cph 
;loadct,0,/silent
;oplot,-reverse(freq_n),ks
;oplot,freq_n,reverse(ks)

;along=interpolate(kw,freq_n*cadence*ntt,ks*pix*nyt)
;along2=interpolate(kw,ntt-freq_n*cadence*ntt,ks*pix*nyt)


;Calculate phase speeds

IF nt mod 2 EQ 0 THEN nt=nt-1 ;deals with even length tracks
icent=mid/2
nlag=nt  ;number of points in cross correlation (make odd)
lag=indgen(nlag)-fix(nlag/2.)
n = mid  ;number of points along track for initial guess (make odd)

num_mc=20 ; number of MC bootstraps to do
randerr=randomn(systime,nt,mid,num_mc) ;generate random errors
lags_pro=fltarr(mid,num_mc,11)
lags_ret=fltarr(mid,num_mc,11)
amplitude2=fltarr(mid,11,2)


;NEEDED FOR RECALCULATING AMPLITUDE
half=(nt-1)/2
X = FINDGEN(half) + 1
IF (nt MOD 2) EQ 0 THEN $
  f = [0.0, X, nt/2, -nt/2 + X]/(nt*cadence) $
ELSE $
  f = [0.0, X, -(nt/2 + 1) + X]/(nt*cadence)

fin=where(f gt 0 and f lt 0.007)
vph=0.390
xi=4.2759
fac=pix/vph/xi


FOR k=5,15 DO BEGIN

    vel=rotate(reform(loop[k,0:nt-1,mid:fpt[1]]),7)
    FOR i=0,mid-1 DO vel[*,i]=vel[*,i]-mean(vel[*,i])

    trans=fft(vel)
    pro_trans=trans             ;select prograde waves (assume nt even, npt odd)
    pro_trans(1:nt/2-1,0:mid/2)=0.
    pro_trans(nt/2+1:nt-1,mid/2+1:mid-1)=0.
    pro_vel=float(fft(pro_trans,/inverse))
    err=fltarr(mid)
    FOR i=0,mid-1 DO err[i]=sqrt(mean((pro_vel[*,i]-smooth(pro_vel[*,i],3,/edge_truncate))^2))
    amplitude2[*,k-5]=sqrt(mean(smooth(pro_vel,[3,1])^2,dim=1))

    errs=randerr*transpose(rebin(err,mid,nt,num_mc),[1,0,2])
    lags_pro[*,*,k-5]=calc_speed(n,mid,nt,icent,pro_vel,errs,num_mc)
    
 
    ret_trans=trans             ;select retrograde waves
    ret_trans(1:nt/2-1,mid/2+1:mid-1)=0.
    ret_trans(nt/2+1:nt-1,0:mid/2)=0.
    ret_vel=float(fft(ret_trans,/inverse))
    FOR i=0,mid-1 DO err[i]=sqrt(mean((ret_vel[*,i]-smooth(ret_vel[*,i],3,/edge_truncate))^2))
    
    errs=randerr*transpose(rebin(err,mid,nt,num_mc),[1,0,2])
    lags_ret[*,*,k-5]=calc_speed(n,mid,nt,icent,ret_vel,errs,num_mc) 

    ;For recalculating amplitude using damping measurements
    amplitude2[*,k-5,0]=sqrt(mean(smooth(pro_vel,[3,1])^2,dim=1))
    FOR j=0,mid-1 DO BEGIN
      dum=abs(fft(smooth(pro_vel[*,j],3)))^2
      amplitude2[j,k-5,1]=sqrt( total( 2.*dum[fin]*exp(2.*fac*j*f[fin]) ) )
    ENDFOR


ENDFOR



fits_p=fltarr(num_mc,5)
fits_r=fltarr(num_mc,5)
resid_p=fltarr(47,num_mc)
resid_r=fltarr(47,num_mc)


FOR i=0,num_mc-1 DO BEGIN

     dist=(findgen(mid)-mid/2)*pix
     mn=(moment(lags_pro[*,i,*],dimension=3))[*,0:1] 
     res=poly_fit(dist,mn[*,0]*cadence,1,measure=sqrt(mn[*,1])*cadence,sigma=sig,chisq=chisq,yfit=fit)
     resid_p[*,i]=mn[*,0]*cadence-fit
     fits_p[i,0:1]=res
     fits_p[i,2:3]=sig
     fits_p[i,4]=chisq

     mn=(moment(lags_ret[*,i,*],dimension=3))[*,0:1] 
     res=poly_fit(dist,mn[*,0]*cadence,1,measure=sqrt(mn[*,1])*cadence,sigma=sig,chisq=chisq,yfit=fit)
     resid_r[*,i]=mn[*,0]*cadence-fit
     fits_r[i,0:1]=res
     fits_r[i,2:3]=sig
     fits_r[i,4]=chisq


ENDFOR

cph_in=1./mean(fits_p[*,1])
cph_in_err=(mean(fits_p[*,3])/mean(fits_p[*,1])^2)
cph_out=1./mean(fits_r[*,1])
cph_out_err=(mean(fits_r[*,3])/mean(fits_r[*,1])^2)

print,'Propagation speed 2nd half', cph_in,cph_out
print,'Propagation speed 2nd half', cph_in_err,cph_out_err
print,sqrt((moment(1./fits_p[*,1]))[1])
print,sqrt((moment(1./fits_r[*,1]))[1])

device,/encapsul,/color,filename=outpath+'lag_example.eps'
loadct,13,/silent
plot,dist,lags_pro[*,0,0]*cadence,xtitle='Distance (Mm)',xst=1,ytitle='Lag (s)'
FOR i=1,10 DO oplot,dist,lags_pro[*,0,i]*cadence
oplot,dist,mean(lags_pro[*,0,*],dimension=3)*cadence,thick=3,color=255
oplot,dist,fits_p[0,0]+fits_p[0,1]*dist,color=150,linestyle=2
device,/close

stop

lags_pro2=fltarr(mid,num_mc,11)
lags_ret2=fltarr(mid,num_mc,11)
amplitude=fltarr(mid,11,3)


FOR k=5,15 DO BEGIN

    vel=reform(loop[k,0:nt-1,fpt[0]:mid-1])
    FOR i=0,mid-1 DO vel[*,i]=vel[*,i]-mean(vel[*,i]) ;remove temporal mean

    trans=fft(vel)
    pro_trans=trans             ;select prograde waves (assume nt even, npt odd)
    pro_trans(1:nt/2-1,0:mid/2)=0.
    pro_trans(nt/2+1:nt-1,mid/2+1:mid-1)=0.
    pro_vel=float(fft(pro_trans,/inverse))
    err=fltarr(mid)
    FOR i=0,mid-1 DO err[i]=sqrt(mean((pro_vel[*,i]-smooth(pro_vel[*,i],3,/edge_truncate))^2))
    


    errs=randerr*transpose(rebin(err,mid,nt,num_mc),[1,0,2])
    lags_pro2[*,*,k-5]=calc_speed(n,mid,nt,icent,pro_vel,errs,num_mc) 
 
    ret_trans=trans             ;select retrograde waves
    ret_trans(1:nt/2-1,mid/2+1:mid-1)=0.
    ret_trans(nt/2+1:nt-1,0:mid/2)=0.
    ret_vel=float(fft(ret_trans,/inverse))
    FOR i=0,mid-1 DO err[i]=sqrt(mean((ret_vel[*,i]-smooth(ret_vel[*,i],3,/edge_truncate))^2))
    
    errs=randerr*transpose(rebin(err,mid,nt,num_mc),[1,0,2])
    lags_ret2[*,*,k-5]=calc_speed(n,mid,nt,icent,ret_vel,errs,num_mc) 

    ;For recalculating amplitude using damping measurements
    amplitude[*,k-5,0]=sqrt(mean(smooth(pro_vel,[3,1])^2,dim=1))
    FOR j=0,mid-1 DO BEGIN
      dum=abs(fft(smooth(pro_vel[*,j],3)))^2
      amplitude[j,k-5,1]=sqrt( total( 2.*dum[fin]*exp(2.*fac*j*f[fin]) ) )
    ENDFOR

ENDFOR


fits_p2=fltarr(num_mc,5)
fits_r2=fltarr(num_mc,5)
resid_p2=fltarr(47,num_mc)
resid_r2=fltarr(47,num_mc)
FOR i=0,num_mc-1 DO BEGIN
 
     dist=(findgen(mid)-mid/2)*pix
     mn=(moment(lags_pro2[*,i,*],dimension=3))[*,0:1] 
     res=poly_fit(dist,mn[*,0]*cadence,1,measure=sqrt(mn[*,1])*cadence,sigma=sig,chisq=chisq,yfit=fit  )
     resid_p2[*,i]=mn[*,0]*cadence-fit
     fits_p2[i,0:1]=res
     fits_p2[i,2:3]=sig
     fits_p2[i,4]=chisq


     mn=(moment(lags_ret2[*,i,*],dimension=3))[*,0:1] 
     res=poly_fit(dist,mn[*,0]*cadence,1,measure=sqrt(mn[*,1])*cadence,sigma=sig,chisq=chisq,yfit=fit)
     resid_r2[*,i]=mn[*,0]*cadence-fit
     fits_r2[i,0:1]=res
     fits_r2[i,2:3]=sig
     fits_r2[i,4]=chisq

ENDFOR
cph_in=1./mean(fits_p2[*,1])
cph_in_err=(mean(fits_p2[*,3])/mean(fits_p2[*,1])^2)
cph_out=1./mean(fits_r2[*,1])
cph_out_err=(mean(fits_r2[*,3])/mean(fits_r2[*,1])^2)

print,'Propagation speed 1st half', cph_in,cph_out
print,'Propagation speed 1st half', cph_in_err,cph_out_err
print,sqrt((moment(1./fits_p2[*,1]))[1])
print,sqrt((moment(1./fits_r2[*,1]))[1])

device,/encapsul,/color,filename=outpath+'prop_speed_hist.eps'
!p.multi=[0,2,2]
loadct,13,/silent
plothist,1./fits_p2[*,1]*1000.,bin=2,xhist,yhist,xtitle='Propagation speed (km s!e-1!n)',charsize=0.85
oplot,xhist,gaussfit(xhist,yhist,nterms=3),linestyle=2,color=255
plots,[cph_in,cph_in]*1000.,[0,max(yhist)]
plothist,1./fits_r2[*,1]*1000.,bin=2,xhist,yhist,xtitle='Propagation speed (km s!e-1!n)',charsize=0.85
oplot,xhist,gaussfit(xhist,yhist,nterms=3),linestyle=2,color=255
plots,[cph_out,cph_out]*1000.,[0,max(yhist)]

plothist,1./fits_p[*,1]*1000.,bin=2,xhist,yhist,xtitle='Propagation speed (km s!e-1!n)',charsize=0.85
oplot,xhist,gaussfit(xhist,yhist,nterms=3),linestyle=2,color=255
plots,[1./mean(fits_p[*,1]),1./mean(fits_p[*,1])]*1000.,[0,max(yhist)]
plothist,1./fits_r[*,1]*1000.,bin=2,xhist,yhist,xtitle='Propagation speed (km s!e-1!n)',charsize=0.85
oplot,xhist,gaussfit(xhist,yhist,nterms=3),linestyle=2,color=255
plots,[1./mean(fits_r[*,1]),1./mean(fits_r[*,1])]*1000.,[0,max(yhist)]

device,/close
!p.multi=0
set_plot,'x'
stop


END