;PURPOSE: - Calculate angular-ly averaged velocity power spectra for CoMP data
;
;
;METHODOLOGY:
;           The aim is to fit the power spectra with a function P=10^af^-b. This will be two
;           parts due to oscillatory bump in power spectra.
;           
;           Using FFT to find periodogram. Following Vaughan, A&A, 431,391 (2005) the routine
;           finds FFT power and converts to log space to make errors homoskedastic.
;           In this case, power in each frequency bin should be close to normal distribution (assuming
;           the sample of the corona we choose all has the same power spectra).
;
;           For each frequency bin, the mean log power and standard error is calculated through 
;           the bootstrap (Enfron) and a Kolmogorov-Smirnoff test is run to check for normality
;           of the distribution of the mean log power. This is done in power_fit_cube. Details of 
;           each frequency distribution are also used.
;
;            
; 
;
;
;INPUTS:- date  - date of aligned cube_ivw
;
;OPTIONAL INPUTS: - apod - degree of apodisation (<1.) - default is 0.1
;                 - extra - string relating to additional characters in cube_ivw name
;                 - /doplot - Outputs diagnostic plots
;                 - phi_ang - select the angle to average over
;                 - rmin - define a minimum radius for averaging
;                 - rmax - define a maximum radius for averaging
;
; EXAMPLE - aa_comp_vel,'20120327',phi_ang=5,rmax=1080,rmin=1050   
;
;
;OUTPUTS: - 
;
;TO DO: - 
;         Problem with 'duff' time-series - current removal of time-series with number of zeros is problematic
;         Might be better to exclude series on rms power values
;


;----------------------------------------------------------------------------
;Calculate r,phi values for image using coord_cart_heilo
function frad_coord,index,image,phi=phi

sz=size(image)                                          
nx=sz(1)
ny=sz(2)      
x=findgen(nx)
x=rebin(x,nx,ny)
y=findgen(ny)   
y=rebin(y,nx,ny)
y=transpose(y)   

coord_cart_helio,index,1.0,x,y,xout,yout,l,b

r=sqrt(xout^2+yout^2)


x_in=(where(xout[*,0]*shift(xout[*,0],1) lt -1))[1]
y_in=(where(xout[*,0]*shift(xout[*,0],1) lt -1))[1]

phi=fltarr(nx,ny)
phi=90-atan(yout/xout)*180./!pi

phi[0:x_in-1,0:y_in-1]=270.-atan(yout[0:x_in-1,0:y_in-1]/xout[0:x_in-1,0:y_in-1])*180./!pi 
phi[0:x_in-1,y_in:ny-1]=270.-atan(yout[0:x_in-1,y_in:ny-1]/xout[0:x_in-1,y_in:ny-1])*180./!pi 

return,r

END

;----------------------------------------------------------------------------
 ; apodizes time-series, with optional de-trending 
FUNCTION apod3dcube,cube,apod,mask,cpg=cpg ;,detrend=detrend
    
  ; get cube dimensions
  sizecube=size(cube)
  nx=sizecube[1]
  ny=sizecube[2]
  nt=sizecube[3]
  apocube = fltarr(nx,ny,nt)

  ; define temporal apodization
  apodt = fltarr(nt)+1
  IF (apod ne 0) THEN BEGIN
    apodrimt = nt*apod
    apodt[0] = (sin(!pi/2.*findgen(apodrimt)/apodrimt))^2
    apodt = apodt*shift(rotate(apodt,2),1)   
  ENDIF
  CPG=total(apodt)/n_elements(apodt) ;Coherent Power Gain - correction factor needed for power after apodisation

  ;temporal mean removal and apodization
  apocube=(cube-rebin(mean(cube,dim=3),nx,ny,nt))*transpose(rebin(apodt,nt,nx,ny),[1,2,0])

  return,apocube
END


;----------------------------------------------------------------------------
; Test for correlations between pixel
; Cross-correlates all pixels and calculates effective number
; of pixels - modified version of Ireland et al., 2015
;
; First version was nearest neighbour:
; Check nearest neighbours
;   0|1|2
;   3| |4
;   5|6|7
;

FUNCTION cor_test,cube,xc,yc,nt,mask,nf 
    
    nelm=n_elements(xc)
    
    rndref=floor(randomu(systime,15000)*nelm)
    ;rnd=floor(randomu(systime,15000)*nelm) ;
    rnd=floor(randomu(systime,15000)*8.)
    ;ccarr=fltarr(nf)

    ;min and max pixel values in angular slice 
    xmm=minmax(xc) 
    ymm=minmax(yc)
    nsmooth=5

    FOR pp=0,14999 DO BEGIN

        ;Check for bad reference pixel
        IF mask[xc[rndref[pp]],yc[rndref[pp]]] EQ 1 THEN BEGIN
            refpix=reform(cube(xc[rndref[pp]],yc[rndref[pp]],*))
;            refpix=reform(cube(xc[rndref[pp]],yc[rndref[pp]],0:nf-1))
           
           ; OLD nearest neighbour version 
            CASE rnd[pp] OF
              0: px=[-1,1]
              1: px=[0,1]
              2: px=[1,1]
              3: px=[-1,0]
              4: px=[1,0]
              5: px=[-1,-1]
              6: px=[0,-1]
              7: px=[1,-1]
            ENDCASE
            xpix=xc[rndref[pp]]+px[0]
            ypix=yc[rndref[pp]]+px[1]
           ;xpix=xc[rnd[pp]]
           ;ypix=yc[rnd[pp]]
           ;ccpix= reform(cube(xpix,ypix,*))
           ccpix=reform(cube[xpix,ypix,0:-1])
        
            ;Check for bad NN pixel and outside range
            gb=mask[xpix,ypix]
            IF GB EQ 1 AND xpix LT xmm[1] AND xpix GT xmm[0] AND ypix LT ymm[1] AND ypix GT ymm[0] THEN BEGIN
                ;Cross-correlate time-series
                tp=max(abs(c_correlate(reform(refpix),reform(ccpix),findgen(nt)-nt/2)))
                IF n_elements(ccarr) eq 0 THEN ccarr=tp ELSE ccarr=[[ccarr],[tp]]                
                
                dis=sqrt((xpix-xc[rndref[pp]])^2+(ypix-yc[rndref[pp]])^2)
                IF n_elements(disarr) eq 0 THEN disarr=dis ELSE disarr=[[disarr],[dis]]  
                ;csp1=smooth(refpix*conj(refpix),nsmooth,/edge_trunc)
                ;csp2=smooth(ccpix*conj(ccpix),nsmooth,/edge_trunc)
                ;cspr=smooth(refpix*conj(ccpix),nsmooth,/edge_trunc)

                ;IF mean(ccarr[*,0]) eq 0 THEN ccarr=real_part(abs(cspr)/sqrt(csp1)/sqrt(csp2)) ELSE $
                 ;                               ccarr=[[ccarr],[real_part(abs(cspr)/sqrt(csp1)/sqrt(csp2))]]
            ENDIF 

        ENDIF 
    ENDFOR

    ;Calculate effective number of pixels
    rho=1.-mean((ccarr))
    neff=1.+rho*(nelm-1)
    
  return,neff
END

FUNCTION comb_plot,meanval,logdat,ix,index,phi_ang,ext_inpath,nt
 
return,!null
END


FUNCTION comb_plot_all,meanval,logdat,ix,index,phi_ang,ext_inpath,freq

IF where(plval eq ix) NE -1 THEN BEGIN
      title=index.date_d$obs+'Angle '+strtrim(string(phi_ang,format='(-F10.2)'),2)+$
      'Frequency '+strtrim(freq[ix],2)+' Hz'

      ;Histogram of mean power
      a=moment(meanval)
      plothist,meanval,xbin,pdf,/noplot
      phisto = PLOT(xbin, pdf, $
      TITLE=title, XTITLE='Log power', YTITLE='Frequency', $
      COLOR='red',/buffer,/stairstep)
      phisto.Refresh,/disable
      t1=text(0.18,0.85,'Mean='+strtrim(tspec[ix,0],2) )
      t2=text(0.18,0.8,'SE='+strtrim(sqrt(tspec[ix,1]),2) )
      t3=text(0.18,0.75,'Skew='+strtrim(a[2],2) )
      t4=text(0.18,0.7,'Kurtosis='+strtrim(a[3],2) )
      t5=text(0.18,0.65,'D!dKS!n='+strtrim(tspec[ix,2],2) )
      t5=text(0.18,0.6,'N!dBS!n='+strtrim(num_bs,2) )
      phisto.Refresh
      phisto.Save, ext_inpath+"/mean_power_plot_"+strtrim(ix,2)+".jpg"
ENDIF

IF where(plval eq ix) NE -1 THEN BEGIN
    title=index.date_d$obs+'Angle '+strtrim(string(phi_ang,format='(-F10.2)'),2)+$
          'Frequency '+strtrim(freq[ix],2)+' Hz'
    ;Histogram of power
    plothist,logdat,xbin,pdf,/noplot
    phisto = PLOT(xbin, pdf, $
        TITLE=title, XTITLE='Log power', YTITLE='Frequency', $
                    COLOR='red',/buffer,/stairstep)
    phisto.Refresh,/disable
    t1=text(0.18,0.85,'Mean='+strtrim(tspec[ix,0],2) )
    t2=text(0.18,0.8,'SD='+strtrim(sqrt(tspec[ix,4]),2) )
    t3=text(0.18,0.75,'Skew='+strtrim(tspec[ix,5],2) )
    t4=text(0.18,0.7,'Kurtosis='+strtrim(tspec[ix,6],2) )
    t5=text(0.18,0.65,'D!dKS!n='+strtrim(tspec[ix,7],2) )
    t5=text(0.18,0.6,'N='+strtrim(tspec[ix,8],2) )
    chk=PLOT(tspec[ix,0]*[1,1],[0,max(pdf)],overplot=phisto,color='red')
    phisto.Refresh
    phisto.Save, ext_inpath+"/power_plot_"+strtrim(ix,2)+".jpg"

    ;Distribution of power as a function of height
    ;Check for obvious trends
    pscatter = PLOT(rval/index.rsun-1., logdat,'b+', xtitle='Solar Radius -1',ytitle='Log Power',$
                                      TITLE='Frequency '+strtrim(freq[ix],2)+' Hz',/buffer)
    pscatter.Save, ext_inpath+"/power_vs_height_"+strtrim(ix,2)+".jpg"
ENDIF  
return,!null
END

;----------------------------------------------------------------------------
;Taken from average_spec.pro and changed
;may cause problems if both are run together
FUNCTION power_fit_cube,nf,input,noise=noise,ext_inpath=ext_inpath, $
                        freq=freq,index=index,rval=rval,doplot=doplot,neff=neff,phi_ang=phi_ang,all=all
    tspec=fltarr(nf,9)
    nelm=n_elements(input[0,*])
    em_corr= 0.57721566 ;Euler constant
    ;which frequency values to plot
    plval=findgen((nf-1)/8)*8

    ;For bootstrap
    num_bs=500
    indx=floor(randomu(systime,nelm,num_bs)*(nelm)) 
    meanval=fltarr(num_bs)
    varval=fltarr(num_bs)

    ;Forget about DC component
    FOR ix=1,nf-1 DO BEGIN

        logdat=alog(reform(input[ix,*])*1d)
        mom=moment(logdat)
        ;upper=mom[0]+3.*sqrt(mom[1])
        ;lower=mom[0]-3.*sqrt(mom[1])
        ;cut=where(logdat gt lower and logdat lt upper )
        ;logdat=logdat[cut]
  
        ;Bootstrap mean distribution
        ;Using KS test to compare to normal CDF
        ;Will need to use Lillefors correction for unknown mean and SD

        FOR bi=0,num_bs-1 DO BEGIN
            meanval[bi]=mean(logdat[indx[*,bi]])
            varval[bi]=(moment(logdat[indx[*,bi]]))[1]
        ENDFOR

        tspec[ix,0]=mean(logdat) + em_corr ;Bias correction (Vaughan 2005)
        tspec[ix,1]=sqrt((moment(meanval))[1]*nelm/neff) ; bootstrap standard error modified by effective pixel
        IF tspec[ix,1] lt !pi/sqrt(6.*nelm) THEN tspec[ix,1]=!pi/sqrt(6.*nelm) 
        testdist=(meanval-tspec[ix,0]+em_corr)/sqrt((moment(meanval))[1])
        ksone,testdist,'gauss_cdf',d,prob
        tspec[ix,2]=d ; test for normality of bootstrap mean
      

        a=moment(logdat)
        tspec[ix,3:6]=a ; Details of PDF of frequency bin
        testdist=(logdat-tspec[ix,0]+em_corr)/sqrt(mean(varval))
        ksone,testdist,'gauss_cdf',d,prob
        tspec[ix,7:8]=[d,nelm] ; test for normality of PDF

        ;Plot mean power
        IF keyword_set(doplot) THEN BEGIN
          nt=n_elements(freq)*2. 
          ;ret=comb_plot(meanval,logdat,ix,index,phi_ang,ext_inpath,nt)
          color=['red','blue','green','cyan','black','orange red','dark blue','purple','yellow']
          title=index.date_d$obs+'Angle '+strtrim(string(phi_ang,format='(-F10.2)'),2)

          IF ix EQ 1 THEN BEGIN
              plothist,meanval,xbin,pdf,/noplot
              denestim = akde(meanval,xbin)
              phisto = PLOT(xbin, denestim, XTITLE='Log power', YTITLE='Frequency', $
              COLOR=color[0],/buffer,xran=[0,-15],title=title)
              phisto.Refresh,/disable

              plothist,logdat,xbin,pdf,/noplot
              phisto_full = PLOT(xbin, pdf, XTITLE='Log power', YTITLE='Frequency', $
              COLOR=color[0],/buffer,/stairstep,xran=[0,-15],title=title)
              phisto_full.Refresh,/disable

              plothist,logdat,xbin,pdf,/noplot
              denestim = akde(logdat,xbin)
              phisto_full_kde = PLOT(xbin, denestim, XTITLE='Log power', YTITLE='Frequency', $
              COLOR=color[0],/buffer,xran=[0,-15],title=title)
              phisto_full_kde.Refresh,/disable
          ENDIF

          IF ix MOD 11 EQ 0 THEN BEGIN
              plothist,meanval,xbin,pdf,/noplot
              denestim = akde(meanval,xbin)
              phisto2 = PLOT(xbin, denestim, XTITLE='Log power', YTITLE='Frequency', $
              COLOR=color[ix/10],overplot=phisto)

              plothist,logdat,xbin,pdf,/noplot
              phisto_full2 = PLOT(xbin, pdf, XTITLE='Log power', YTITLE='Frequency', $
              COLOR=color[ix/10],overplot=phisto_full,/stairstep)

              plothist,logdat,xbin,pdf,/noplot
              denestim = akde(logdat,xbin)
              phisto_full_kde2 = PLOT(xbin, denestim, XTITLE='Log power', YTITLE='Frequency', $
              COLOR=color[ix/10],overplot=phisto_full_kde)
          ENDIF
          IF ix EQ nt/2-1 THEN BEGIN
            phisto.Refresh
            phisto.Save, ext_inpath+"/mean_power_plot.jpg"
            phisto_full.Refresh
            phisto_full.Save, ext_inpath+"/power_plot.jpg"
            phisto_full_kde.Refresh
            phisto_full_kde.Save, ext_inpath+"/power_plot_kde.jpg"
          ENDIF  
          IF keyword_set(all) THEN BEGIN
              ret=comb_plot_all(meanval,logdat,ix,index,phi_ang,ext_inpath,freq)
          ENDIF 
        ENDIF

      ;stop
      ENDFOR
  
    return,tspec

END

;----------------------------------------------------------------------------
FUNCTION fit_spec,instr,nf,cut_off=cut_off,doplot=doplot,ext_inpath=ext_inpath,all=all,postfif=postfif

spectra=instr.spectra[*,0]
freq=instr.freq
errors=instr.spectra[*,1]  ;Assumes Standard error on mean is Gaussian distributed
date=instr.date

x=alog(freq[1:nf-2]) ; Avoid fitting of Nyquist frequency - different uncertainty to
                     ; to other frequency bins (Vaughan 2005)


pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},3)
pi[0].limited(0)=1 & pi[0].limits(0)=1e-8
pi[1].limited(1)=1 & pi[1].limits(1)=-0.1
pi[2].limited(0)=1 & pi[2].limits(0)=5e-6

start=double([2e-6,-1.2,5e-3])

resp=mpfitfun('mypowline',x,spectra[1:nf-2,0],errors[1:nf-2],start,/quiet,$
perror=perror,bestnorm=bestnorm,dof=dof,parinfo=pi,nfree=nfree)
instr.fit_line=resp
instr.err_line=perror
instr.chisq_line=bestnorm/dof
;Calculate log likelihood and Akaike Information Criterion
instr.loglik_line=-(nf-2)*alog(2.*!pi)-total(alog(errors[1:nf-2]))-bestnorm/2.
instr.aic_line=2.*nfree-2*instr.loglik_line+2.*nfree*(nfree-1)/(dof-1)

;Set limits for fitting - required to stop NANs during function evaluation
;occurs due to alog in mypowgauss
pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},6)
pi[0].limited(0)=1 & pi[0].limits(0)=1e-8
pi[1].limited(1)=1 & pi[1].limits(1)=-0.1
pi[2].limited(0)=1 & pi[2].limits(0)=5e-6
pi[3].limited(0)=1 & pi[3].limits(0)=1e-5
pi[4].limited([0,1])=1 & pi[4].limits(0)=-7. & pi[4].limits(1)=-5 
pi[5].limited(1)=1 & pi[5].limits(1)=0.8

IF keyword_set(postfif) THEN startvals=double([2e-5,-1.2,5e-3,0.01,-5.5,0.4]) $
                        ELSE startvals=double([2e-6,-1.2,5e-3,0.001,-5.5,0.4])

respg=mpfitfun('mypowgauss',x,spectra[1:nf-2,0],errors[1:nf-2],startvals,$
perror=perror,bestnorm=bestnorm,dof=dof,parinfo=pi,/quiet,nfree=nfree)
instr.fit_pg=respg
instr.err_pg=perror
instr.chisq_pg=bestnorm/dof
;Calculate log likelihood and Akaike Information Criterion
instr.loglik_pg=-(nf-2)*alog(2.*!pi)-total(alog(errors[1:nf-2]))-bestnorm/2.
instr.aic_pg=2.*nfree-2*instr.loglik_pg+2.*nfree*(nfree-1)/(dof-1)

        ;Plot fit results
        IF keyword_set(doplot) THEN BEGIN
           title=instr.date+' Angle '+strtrim(string(instr.phi,format='(-F10.2)'),2)
           p = ERRORPLOT(x, spectra[1:nf-2,0],errors[1:nf-2], $
                      TITLE=title, YTITLE='Log!de!n power', XTITLE='Log!de!n Frequency', $
                      layout=[1,2,1],position=[0.1,0.4,0.9,0.9], /buffer, NAME='Fourier spectra',line=6,$
                      symbol='+')
           p.Refresh,/disable

           model=mypowline(x,resp)
           modplot=plot(x,model,overplot=p,COLOR='blue', NAME='Model 1',thick=2)
           resid2=(spectra[1:nf-2,0]-model)/errors[1:nf-2]

           model=mypowgauss(x,respg)
           modplot2=plot(x,model,overplot=p,COLOR='red', NAME='Model 2',thick=2)
           resid=(spectra[1:nf-2,0]-model)/errors[1:nf-2]

           
           
           residplot=PLOT(x,resid,'r',layout=[1,2,2],position=[0.1,0.1,0.9,0.3],/current,xtitle='Log!de!n Frequency',$
                  ytitle='Normalised residuals')
           residplot2=PLOT(x,resid2,'b--',overplot=residplot)

           t1=text(0.67,0.84,'$\chi^2_r$='+strtrim(string(instr.chisq_line,format='(-F10.2)'),2) +' Model 1')
           t2=text(0.67,0.79,'$\chi^2_r$='+strtrim(string(instr.chisq_pg,format='(-F10.2)'),2) +' Model 2' )
           t3=text(0.67,0.74,'$\Delta_{AIC}$='+strtrim(fix(instr.aic_line-instr.aic_pg),2) )
           line=plot([-9,-4],[0,0],'--',trans=50,overplot=residplot)

           leg = LEGEND(TARGET=[p,modplot,modplot2], POSITION=[0.38,0.54], $
                      /AUTO_TEXT_COLOR, font_size=8)
           p.Refresh
           p.Save, ext_inpath+"/log_power_fit_"+strtrim(instr.phi,2)+".jpg"
          
           plothist,resid,xhist,pdf,/noplot,bin=0.5
           normplot=plot(xhist,pdf,'r',/stairstep,title='Normalised residuals',/buffer)
           normplot.Refresh,/disable
          
           plothist,resid2,xhist,pdf,/noplot,bin=0.5
           normplot2=plot(xhist,pdf,'b--',/stairstep,title='Normalised residuals',overplot=normplot)
           gx=findgen(20)/2-5.
           gaussplot=plot(gx,max(pdf)*exp(-gx^2/2.),'k-:',overplot=normplot)
           normplot.Refresh
           normplot.Save, ext_inpath+"/log_power_fit_NR_"+strtrim(instr.phi,2)+".jpg"
           
        ENDIF

return,!null

END


;----------------------------------------------------------------------------
PRO AA_comp_vel,date,index=index,$           
                apod=apod,extra=extra,doplot=doplot,$
                do_atrous=do_atrous,phi_ang=phi_ang,rmin=rmin,rmax=rmax,all=all, $
                postfif=postfif


; paths where the CoMP data is located and
; where the .sav-files should be written to
inpath  = 'analysis/CoMP/wave_tracking_output/'+date+'/'
outpath = 'analysis/CoMP/wave_tracking_output/'+date+'/'


resolve_routine, 'gauss_cdf',/is_function
resolve_routine, 'kde',/is_function

IF NOT keyword_set(extra) THEN restore,inpath+'cube_ivw_'+date+'.sav',/verbose $
ELSE restore,inpath+'cube_ivw_'+date+'_'+extra+'.sav',/verbose

data=cube_v
sz=size(data)
nx=sz(1) & ny=sz(2) & nt=sz(3)
xscale=index.cdelt1 ;pixel size in arcsec

;Make series length even to keep things easy!!
IF nt mod 2 EQ 1 THEN BEGIN 
     nt=nt-1
     data=temporary(data[0:-1,0:-1,0:nt-1])
ENDIF

; Upper and lower limit for the radius in
; !be aware that the maximum offset of the cross-correlation
; (if performed) is added to this value!
IF n_elements(rmin) EQ 0 THEN lower_r = index.lower_r*xscale ELSE lower_r=rmin
IF n_elements(rmax) EQ 0 THEN upper_r = index.upper_r*xscale ELSE upper_r=rmax
maxoffset=index.maxoffset*xscale

; Create arrays of r & phi values
r=frad_coord(index,data[*,*,0],phi=phi)

; create mask & mask out pixels not to invert
mask=index.mask
mask[where(r lt lower_r)]=0
mask[where(r gt upper_r)]=0
data_mask=data
for i=0,nt-1 do data_mask[*,*,i]=data[*,*,i]*mask

;Apodise the time-series for each pixel
print,'Starting apodisation'
IF n_elements(apod) LT 1 THEN apod=0.1
apocube=apod3dcube(data_mask,apod,mask,cpg=cpg)


;Very high pass filter
;Good at isolating noise and keeping small scales features
IF keyword_set(do_atrous) THEN BEGIN
   FOR i=0,nt-1 DO BEGIN
       atrous,apocube[*,*,i],decomposition=a,n_scales=1
       oas=a[*,*,1]-smooth(a[*,*,1],3)
       apocube[*,*,i]=temporary(apocube[*,*,i])-oas
   ENDFOR
ENDIF

;Estimates the error levels for power per frequency bin
;using Parcivel's theorem
pwr_err_estim=total((apocube-smooth(apocube,[1,1,3],/edge_truncate))^2,3)/nt^2

;RMS value
v_err_rms=sqrt(mean( (apocube-smooth(apocube,[1,1,3],/edge_truncate))^2,dim=3) )

print,'Starting Fourier transform'
dft=fft(apocube,dim=3)


X = FINDGEN((nt - 1)/2) + 1
freq = [0.0, X, nt/2]/nt/index.norm_cadence 
nf=n_elements(freq)

;Old code - should be removed
cut_off=where(freq gt 1./(3.*index.norm_cadence)) ; frequency elements less than 90s

;Update mask to eliminate pixels from analysis with either excessive or minuscule power values
;Based on full FOV values
mntemp=moment(alog(pwr_err_estim[where(pwr_err_estim ne 0)]))
lower_cut=exp(mntemp[0]-5.*sqrt(mntemp[1]))
upper_cut=exp(mntemp[0]+5.*sqrt(mntemp[1]))
mask[where(pwr_err_estim lt lower_cut)]=0
mask[where(pwr_err_estim gt upper_cut)]=0
for i=0,nt-1 do apocube[*,*,i]=apocube[*,*,i]*mask

;Phi averaging
IF n_elements(phi_ang) EQ 0 THEN phi_ang=45
no_phi=360/phi_ang

;data store
temp={date:index.date_d$obs,freq:fltarr(nf),spectra:fltarr(nf,9),phi:0.0,sub:fltarr(cut_off[0]),$
            fit_pg:fltarr(6),err_pg:fltarr(6),chisq_pg:0.0,loglik_pg:0.0,aic_pg:0.0,fit_line:fltarr(3),$
            err_line:fltarr(3),chisq_line:0.0,loglik_line:0.0,aic_line:0.0}

temp.freq=freq

;----------------------------------------------
;Plotting of image and angle contours
;Define colour table
aia_lct,wavelnth='193',/load
tvlct,red,green,blue,/get
cols=[[red],[green],[blue]]

yran=([0,ny-1]-index.crpix2)*xscale
xran=([0,nx-1]-index.crpix1)*xscale
intplot=image(reform(cube_i[*,*,0]>0<25.),axis_style=0,rgb_table=cols,title=index.date_d$obs)
intplot.Refresh,/disable
xaxis=axis('x',location=0,coord_trans=[xran[0],xscale],title='Distance (arcsec)')
xaxis=axis('y',location=0,coord_trans=[yran[0],xscale],title='Distance (arcsec)')

cvalue=phi_ang*findgen(no_phi)
contplot=contour(phi,overplot=intplot,c_value=cvalue[0:-1:2],c_linestyle=2,color='white',c_label_show=0)
contplot2=contour(mask,overplot=intplot,c_value=[1],c_linestyle=2,color='white',c_label_show=0)

xloc=index.crpix1+(280.)*sin(cvalue[0:-1:2]*!const.dtor)
yloc=index.crpix2+280.*cos(cvalue[0:-1:2]*!const.dtor) 
mylabels=text(xloc,yloc,strtrim(uint(cvalue[0:-1:2]),2),overplot=contplot,color='white',/data)
contplot.Scale,0.9,0.9
intplot.Refresh
contplot.Save,inpath+'angle_plot.eps'
;---------------------------------------------

print,'Do you want to continue? (.c)'
stop

print,'Starting angle averaging'
FOR i=1,no_phi-1 DO BEGIN
     phold=!null
     rval=!null

     ;Determine pixels that lie between two angles          
     wedge=where(phi gt phi_ang*(i-1) and phi lt phi_ang*i and mask gt 0)
     xc=wedge mod nx
     yc=wedge/nx
     
     
     FOR j=0,n_elements(xc)-1 DO BEGIN
         IF mask[xc(j),yc(j)] eq 1. THEN BEGIN  
            
            dft_temp=reform(dft[xc(j),yc(j),0:nf-1]) 
            
            IF n_elements(phold) EQ 0 THEN phold=2.*abs(dft_temp)^2/cpg^2 ELSE $
            phold=[[phold],[2.*abs(dft_temp)^2/cpg^2]]

            IF n_elements(rval) EQ 0 THEN rval=r[xc[j],yc[j]] ELSE $
                                       rval=[[rval],[r[xc[j],yc[j]] ]]
            IF n_elements(rerr) EQ 0 THEN rerr=pwr_err_estim[xc[j],yc[j]] ELSE $
                                        rerr=[[rerr],[pwr_err_estim[ xc[j],yc[j]] ]]
         ENDIF 
     ENDFOR
     
     
    ;Make sure there are more than 10 time-series before running     
     IF n_elements(phold) GT nf*10 THEN BEGIN
     ;correlation test following Ireland et al (2015)
     neff=cor_test(apocube,xc,yc,nt,mask,nf)

            ;Creates folders for plotting
            IF keyword_set(doplot) THEN BEGIN
                IF strlen(file_search(inpath+'aa_diagnostics')) EQ 0 THEN spawn,'mkdir '+inpath+ 'aa_diagnostics'
       
                ext_inpath=inpath+'aa_diagnostics/'+strtrim(phi_ang*(i-1),2)
                IF strlen(file_search(ext_inpath)) EQ 0 THEN $
                                spawn,'mkdir '+inpath+ 'aa_diagnostics/'+strtrim(phi_ang*(i-1),2)
            ENDIF


            ;Finds mean values for each frequency bin
            res=power_fit_cube(nf,phold,noise=noise,ext_inpath=ext_inpath, $
                            freq=freq,index=index,rval=rval,neff=neff,phi_ang=phi_ang*(i-1),doplot=doplot,all=all)
           
            temp.spectra=res
            temp.phi=phi_ang*(i-1)

            ;Fit mean values
            res=fit_spec(temp,nf,cut_off=cut_off,doplot=doplot,ext_inpath=ext_inpath,postfif=postfif)

           IF n_elements(results) EQ 0 THEN results=temp ELSE results=[temporary(results),temp]
    
      ENDIF 
      counter,i,no_phi, 'Calculating spectra '

ENDFOR
save,results,filename=outpath+'aa_vel_results.idlsav'



END