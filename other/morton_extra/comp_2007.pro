
;For radial average of CoMP 2007 data power spectra
;Nat pub 2017
PRO comp_2007,reb=reb



restore,'/Users/richardmorton/analysis/comp/2007/velocity_interp.sav'


velocity=velocity[*,*,0:350]
sz=size(velocity)
nx=sz(1)
ny=sz(2)
nt=sz(3)

IF keyword_set(reb) THEN velocity=rebin(velocity,nx/reb,ny/reb,nt) ELSE reb=1
IF keyword_set(reb) THEN extra='rebin_'+strtrim(reb,2) ELSE extra=''
sz=size(velocity)
nx=sz(1)
ny=sz(2)

low=240/reb
high=260/reb
print,'Spatial sampling (solar rad)',4.46*reb*725./695e3



tvim,velocity[*,*,0]
;print,'pick points on occulting disk edge'

pick,x1,y1

pick,x2,y2


;Calculate center circle
qc=sqrt((x2-x1)^2+(y2-y1)^2)

x3=min([x1,x2])+abs(x2-x1)/2. & y3=min([y1,y2])+abs(y2-y1)/2.

xc = x3 - sqrt((228./reb)^2-(qc/2)^2)*(y1-y2)/qc
yc = y3 - sqrt((228./reb)^2-(qc/2)^2)*(x2-x1)/qc

IF xc lt 0 THEN BEGIN
xc = x3 + sqrt((228./reb)^2-(qc/2)^2)*(y1-y2)/qc
yc = y3 + sqrt((228./reb)^2-(qc/2)^2)*(x2-x1)/qc
ENDIF



xarr=rebin(findgen(nx),nx,ny)
yarr=fltarr(nx,ny)
for i=0,nx-1 do yarr[i,*]=findgen(ny)
r=sqrt((xarr-xc)^2+(yarr-yc)^2)

contour,r,levels=low,/over
contour,r,levels=high,/over

in=where(r gt low and r lt high)
print,'Number of elements',n_elements(in)

series=fltarr(n_elements(in),nt)

apocube = fltarr(nx,ny,nt)
apod=0.1

; define temporal apodization
apodt = fltarr(nt)+1
IF (apod ne 0) THEN BEGIN
    apodrimt = nt*apod
    apodt[0] = (sin(!pi/2.*findgen(apodrimt)/apodrimt))^2
    apodt = apodt*shift(rotate(apodt,2),1)   
ENDIF
CPG=total(apodt)/n_elements(apodt) ;Coherent Power Gain - correction factor needed for power after apodisation

 ;temporal mean removal and apodization
 apocube=(velocity-rebin(mean(velocity,dim=3),nx,ny,nt))*transpose(rebin(apodt,nt,nx,ny),[1,2,0])

 print,'Doing Fourier'
 ftcube=2.*abs(fft(apocube,dim=3))^2/cpg^2

FOR i=0,nt-1 DO BEGIN 
    temp=ftcube[*,*,i]
	series[*,i]=temp[in]
ENDFOR

X = FINDGEN((nt - 1)/2) + 1
freq = [0.0, X, nt/2]/nt/29.
nf=n_elements(freq)
ave_spec=fltarr(nf,2)



nelm=n_elements(series[*,0])
neff=nelm
em_corr= 0.57721566 ;Euler const

;For bootstrap
num_bs=500
indx=floor(randomu(systime,nelm,num_bs)*(nelm)) 
meanval=fltarr(num_bs)
varval=fltarr(num_bs)

print,'Calculating average'
;Forget about DC component
FOR ix=1,nf-1 DO BEGIN
     logdat=alog(reform(series[*,ix])*1d)
      mom=moment(logdat)

      ;Bootstrap mean distribution
      ;Using KS test to compare to normal CDF
      ;Will need to use Lillefors correction for unknown mean and SD
       FOR bi=0,num_bs-1 DO BEGIN
           meanval[bi]=mean(logdat[indx[*,bi]])
           varval[bi]=(moment(logdat[indx[*,bi]]))[1]
       ENDFOR

      ave_spec[ix,0]=mean(logdat) + em_corr ;Bias correction (Vaughan 2005)
      ave_spec[ix,1]=sqrt((moment(meanval))[1]*nelm/neff) ; bootstrap standard error modified by effective pixel
       IF ave_spec[ix,1] lt !pi/sqrt(6.*nelm) THEN ave_spec[ix,1]=!pi/sqrt(6.*nelm) 
ENDFOR

spec2007={freq:freq, spec:ave_spec}
save,spec2007,filename='analysis/comp/2007/coronal_average_spec'+extra+'.idlsav'
stop 
END