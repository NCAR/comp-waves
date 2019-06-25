;PURPOSE: Calculates the errors for the analytic fits to the line profile.
;         Analytic profile is given in Tian et al., 2012
;
;INPUT: date - date of observations
;OPTIONAL INPUTS: see - set to include a value of seeing uncertainties
;                       e.g. see=2e-3 - the value has to be the variance (sigma_see^2)
;
;OUTPUTS: dv - errors on Doppler velocity
;         dw - errors on Doppler width
;         di - errors on central intensity
;         Note - if error values eq 0 THEN fit is bad
;
;
;HISTORY: Ver 1. Created RJM 14/07/2014
;         Ver 2. Seperated out Gaussian & Poisson noise
;
;TO DO: Add background intensity errors + read noise + dark + flat
;
;
;RESTRICTIONS: Assumes all data in sequence is averaged over the same number of frames
;              i.e. goes off frame 0



PRO error_on_fits_v2,date,dv=dv,dw=dw,di=di,see=see


; paths where the CoMP level 1 data is located
inpath = 'data/CoMP/'+strmid(date,8,4,/reverse)+'/'+strmid(date,3,2,/reverse)+'/'+strmid(date,1,2,/reverse)+'/level1'

outpath = 'data/CoMP/wave_tracking_output/'+date+'/'


files=find_files('*.fts.gz',inpath)

;SCALING FOR PHOTONS
pscale=875.

;Read in headers for initial configuration
index1=headfits(files[0],ext=1)
index2=headfits(files[0],ext=2)
index3=headfits(files[0],ext=3)





;Get wavelength separation
lam1=sxpar(index1,'waveleng')
lam2=sxpar(index2,'waveleng')
del=(lam2-lam1)/lam2*2.99e8/1000. ;d in km/s


;Get number of frames averaged together
frame_av=sxpar(index2,'naverage')*1d

;Get array dimensions
nx=sxpar(index1,'naxis1')
ny=sxpar(index1,'naxis2')
nt=160; n_elements(files)

; create mask
; mask out pixels to invert
;lower_r = 228.
;upper_r = 288.
;mask=intarr(nx,ny)
;x=rebin(findgen(ny)-(nx/2.)+0.5,nx,ny)
;y=transpose(x)
;r=sqrt(x^2+y^2)
;good=where(r ge lower_r and r le upper_r)
;mask(good)=1.



dv=fltarr(nx,ny,nt,2)
dw=fltarr(nx,ny,nt,2)
di=fltarr(nx,ny,nt,2)

delI1=fltarr(nx,ny,2)
delI2=fltarr(nx,ny,2)
delI3=fltarr(nx,ny,2)

FOR i=0,nt-1 DO BEGIN
    
    ;Load intensity at each wavelength position
    I1=readfits(files(i),ext=1,/silent)
    I2=readfits(files(i),ext=2,/silent)
    I3=readfits(files(i),ext=3,/silent)

    ;Calculate Doppler velocity and widths
    analytic_gauss_fit_fast,I1,I2,I3,del,v,w,ic,/get_rid
    ic=ic*875.

   ;Get rid of values lt 0
    I1=pscale*temporary(I1)>(1e-4)
    I2=pscale*temporary(I2)>(1e-4)
    I3=pscale*temporary(I3)>(1e-4)

    ;Background images for each wavelength position
    bck1=pscale*readfits(files(i),ext=10,/silent)
    bck2=pscale*readfits(files(i),ext=11,/silent)
    bck3=pscale*readfits(files(i),ext=12,/silent)

    

    ;Readnoise
    sig_r=70.

   
    ;Error on intensity values
    ;Separate Gaussian & poisson values for 
    ;noise simulations 
    IF n_elements(see) GT 0 THEN BEGIN
       delI1[*,*,0]=sqrt(sig_r^2+sobel(I1)^2*see)/sqrt(frame_av)
       delI1[*,*,1]=sqrt(I1+bck1)/sqrt(frame_av)
       delI2[*,*,0]=sqrt(sig_r^2+sobel(I2)^2*see)/sqrt(frame_av)
       delI2[*,*,1]=sqrt(I2+bck2)/sqrt(frame_av)
       delI3[*,*,0]=sqrt(sig_r^2+sobel(I3)^2*see)/sqrt(frame_av)
       delI3[*,*,1]=sqrt(I3+bck3)/sqrt(frame_av)
    ENDIF ELSE BEGIN
      delI1=sqrt(I1+bck1+sig_r^2)/sqrt(frame_av)
      delI2=sqrt(I2+bck2+sig_r^2)/sqrt(frame_av)
      delI3=sqrt(I3+bck3+sig_r^2)/sqrt(frame_av)
    ENDELSE
    

    ;##################################
    ;ANALYTIC FORMULA FOR GAUSSIAN FITS
    ;##################################
    ;v=-0.5*del*(a-b)/(a+b)
    ;w=(-2*del^2/(a+b))^0.5
    ;I0=I2*exp(v^2/w^2)
    ;
    ;a=ln(I3/I2)
    ;b=ln(I1/I2)
    ;del=lamda2-lambda1 (lambda3-lambda2)

    ;##################################
    ;Error on Doppler shift
    ;v=-0.5*d*[ln(I3)-ln(I1)]/(ln(I3)+ln(I1)-2ln(I2))
    ;dv^2=sum(dv/dI_j*delI_j)^2 ; Standard error propagation
    ;
    ;##################################

    sdel=(alog(I3)+alog(I1)-2.*alog(I2))

    dvdI1=0.5/I1/sdel + 0.5*(alog(I3)-alog(I1))/I1/sdel^2
    dvdI2=-(alog(I3)-alog(I1))/I2/sdel^2
    dvdI3=-0.5/I3/sdel + 0.5*(alog(I3)-alog(I1))/I3/sdel^2

    dvsq_g=(del*dvdI1*delI1[*,*,0])^2+(del*dvdI2*delI2[*,*,0])^2+(del*dvdI3*delI3[*,*,0])^2
    dvsq_p=(del*dvdI1*delI1[*,*,1])^2+(del*dvdI2*delI2[*,*,1])^2+(del*dvdI3*delI3[*,*,1])^2
    ;Get rid of -NaN's
    dvsq_g[where(finite(dvsq_g,/nan,sign=-1))]=0.
    dvsq_p[where(finite(dvsq_p,/nan,sign=-1))]=0.
    
    dv[*,*,i,0]=sqrt(dvsq_g)
    dv[*,*,i,1]=sqrt(dvsq_p)
    

    ;##################################
    ;Error on Doppler width
    ;w=[2*d^2/(-ln(I3)-ln(I1)+2ln(I2))]^0.5
    ;dw^2=sum(dw/dI_j*delI_j)^2 ; Standard error propagation
    ;
    ;##################################
    sdel=(-alog(I3)-alog(I1)+2*alog(I2))

    dwdI1=0.5*sqrt(2*del^2)/I1/sdel^(1.5)
    dwdI2=-sqrt(2*del^2)/I2/sdel^(1.5)
    dwdI3=0.5*sqrt(2*del^2)/I3/sdel^(1.5)

    dwsq_g=(dwdI1*delI1[*,*,0])^2+(dwdI2*delI2[*,*,0])^2+(dwdI3*delI3[*,*,0])^2
    dwsq_p=(dwdI1*delI1[*,*,1])^2+(dwdI2*delI2[*,*,1])^2+(dwdI3*delI3[*,*,1])^2
    
    ;Get rid of -NaN's
    dwsq_g[where(finite(dwsq_g,/nan,sign=-1))]=0.
    dwsq_p[where(finite(dwsq_p,/nan,sign=-1))]=0.
    
    dw[*,*,i,0]=sqrt(dwsq_g)
    dw[*,*,i,1]=sqrt(dwsq_p)

    ;##################################
    ;Error on central intensity
    ;i=1./(8*del)*(lnI3-lnI1)^2/(ln(I3)+ln(I1)-2ln(I2))
    ;di^2=sum(di/dI_j*delI_j)^2 ; Standard error propagation
    ;
    ;##################################
    sdel=(alog(I3)+alog(I1)-2*alog(I2))

    didI1=1./4.*(alog(I3)-alog(I1))/I1/sdel+1./8.*(alog(I3)-alog(I1))^2/I1/sdel^2
    didI2=1.-1./4.*(alog(I3)-alog(I1))^2/I2/sdel^2
    didI3=-1./4.*(alog(I3)-alog(I1))/I3/sdel+1./8.*(alog(I3)-alog(I1))^2/I3/sdel^2

    disq_g=(ic*didI1*delI1[*,*,0])^2+(exp(v^2/w^2)*didI2*delI2[*,*,0])^2+(ic*didI3*delI3[*,*,0])^2
    disq_p=(ic*didI1*delI1[*,*,1])^2+(exp(v^2/w^2)*didI2*delI2[*,*,1])^2+(ic*didI3*delI3[*,*,1])^2

    ;Get rid of -NaN's
    disq_g[where(finite(disq_g,/nan,sign=-1))]=0.
    disq_p[where(finite(disq_p,/nan,sign=-1))]=0.

    di[*,*,i,0]=sqrt(disq_g)
    di[*,*,i,1]=sqrt(disq_p)

    counter,i,nt,'Stuff',/percent

ENDFOR

IF n_elements(see) GT 0 THEN save,dv,dw,di,filename=outpath+'/cube_ivw_errors_incsee.sav' $
ELSE save,dv,dw,di,filename=outpath+'/cube_ivw_errors.sav'

;!p.multi=[0,3,1]
;in=where(ic gt 1. and ic lt 40.*875.)
;plothist,(di[in]/ic[in])<0.2,xrange=[0.0,0.2],bin=0.001
;plothist,(dv[in]/v[in]) <1,xrange=[0.0,1],bin=0.005
;plothist,(dw[in]/w[in]) <0.3,xrange=[0.0,0.3],bin=0.005
;!p.multi=0

END