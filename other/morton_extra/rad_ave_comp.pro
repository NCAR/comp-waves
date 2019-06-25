;PURPOSE: - Calculate radially and angular-ly averaged power spectra for CoMP data
;
;INPUTS:- index - initial of observations (only used to calculate the radial coordinates) 
;         date  - date of aligned cube_ivw
;
;OPTIONAL INPUTS: - maxoffset - maximum offset from cross-correlation in wave_tracking.pro
;                 - apod - degree of apodisation (<1.) - default is 0.1
;
;                   extra - string relating to additional characters in cube_ivw name
;                   set - use to specify which CoMP observable to use, default is velocity
;                         use set='dw' or 'di' to use Doppler width or Intensity
;                   /errors - will find error file for the day and plot power spectra for the errors
;                   /do_phi - will provide the data averaged over different phi angles (no radial seperation)
;
;
;OUTPUTS: - ravpow - radially averaged power array of frequency vs height
;           phipow - angular-ly averaged power array of frequency vs angle
;           rplot  - array of radial coordinates in arcsec
;
;



;Calcuate r,phi values for image using coord_cart_heilo
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
function linear, x, p            ;used in temporal detrending 
  ymod = p[0] + x * p[1]
  return, ymod
end


; apodizes time-series, with optional de-trending 
;----------------------------------------------------------------------------
FUNCTION apod3dcube,cube,apod ;,detrend=detrend
    
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

  ; temporal detrending, not per px, only mean-image trend 
  ttrend = fltarr(nt)
  tf = findgen(nt) + 1
  FOR ix=0, nx-1 DO BEGIN
    FOR iy=0, ny-1 DO BEGIN
        ts = cube[ix,iy,0:nt-1]
        res=poly_fit(findgen(nt),ts,1,yfit=fit)
        apocube[ix,iy,*]=(ts-fit)*apodt
        
    ENDFOR
  ENDFOR
 
  ; done
  return,apocube
END






;----------------------------------------------------------------------------
PRO rad_ave_comp,date,index=index,ravpow=ravpow,phiavpow=phiavpow,maxoffset=maxoffset,$           
                            rplot=rplot,apod=apod,extra=extra,set=set,errors=errors,$
                            do_phi=do_phi,do_atrous=do_atrous


; paths where the CoMP data is located and
; where the .sav-files should be written to
inpath  = 'data/CoMP/wave_tracking_output/'+date+'/'
outpath = 'data/CoMP/wave_tracking_output/'+date+'/'

inpath_ind  = 'data/CoMP/'+strmid(date,8,4,/reverse)+'/'+strmid(date,3,2,/reverse)+'/'+strmid(date,1,2,/reverse)+'/'
files=find_files('*.fts.gz',inpath_ind)
mreadfits,files(0),index


IF keyword_set(errors) THEN restore,inpath+'cube_ivw_errors.sav',/verbose ELSE $
IF NOT keyword_set(extra) THEN restore,inpath+'cube_ivw_'+date+'.sav',/verbose $
ELSE restore,inpath+'cube_ivw_'+date+'_'+extra+'.sav',/verbose


IF n_elements(set) EQ 0 THEN set ='dv'

IF keyword_set(errors) THEN BEGIN
    CASE SET OF
       'dv':data=dv<100.
       'dw':data=dw
       'di':data=di
    ENDCASE
ENDIF ELSE BEGIN 
    CASE SET OF
       'dv':data=cube_v
       'dw':data=cube_w
       'di':data=cube_i
    ENDCASE
ENDELSE


sz=size(data)
nx=sz(1)
ny=sz(2)
nt=sz(3)
xscale=index.cdelt1 ;pixel size in arcsec


; Upper and lower limit for the radius in
; !be aware that the maximum offset of the cross-correlation
; (if performed) is added to this value!
lower_r = 228.*xscale
upper_r = 288.*xscale

IF n_elements(maxoffset) LT 1 THEN maxoffset=0. ELSE maxoffset=maxoffset*xscale
IF n_elements(apod) LT 1 THEN apod=0.1


image=data[*,*,0]

r=frad_coord(index,image,phi=phi)

; create mask
; mask out pixels to invert
mask=intarr(nx,ny)

good=where(r ge lower_r+maxoffset and r le upper_r)
mask(good)=1.

data_mask=data
for i=0,nt-1 do data_mask[*,*,i]=data[*,*,i]*mask



no_r=floor((upper_r-lower_r-maxoffset)/xscale)
rplot=findgen(no_r)*xscale+lower_r+maxoffset

ravpow=fltarr(nt/2,no_r,2)
rpow_elem=0.

print,'Starting'

;Apodise the time-series for each pixel
apocube=apod3dcube(data_mask,apod)


;Very high pass filter
;Good at isolating noise and keeping small scales features
;
IF keyword_set(do_atrous) THEN BEGIN
   FOR i=0,nt-1 DO BEGIN
       atrous,apocube[*,*,i],decomposition=a,n_scales=1
       oas=a[*,*,1]-smooth(a[*,*,1],3)
       ;oas=oas-smooth(oas,3)
       apocube[*,*,i]=temporary(apocube[*,*,i])-oas
   ENDFOR
ENDIF


;Radial averaging
FOR i=1,no_r DO BEGIN

   ring=where(r gt (lower_r+maxoffset+(i-1)*xscale) and r lt (lower_r+maxoffset+i*xscale) )
 ;print,lower_r+maxoffset+(i-1)*xscale,lower_r+maxoffset+(i)*xscale
 ;print,n_elements(ring)
  xc=ring mod nx
  yc=ring/nx

  phold=fltarr(nt/2)
  FOR j=0,n_elements(xc)-1 DO BEGIN
  
     no_zero=n_elements(where(apocube[xc(j),yc(j),*] eq 0.))
      
     IF no_zero lt 5 THEN BEGIN 
        times=apocube[xc(j),yc(j),*]
        

        ;Excludes pixels with time-series that have spikes greater
        ;than pm 4 sigma from the mean
        
        vals=moment(times)
        cl=(vals[0]+4.*sqrt(vals[1]))
        cll=(vals[0]-4.*sqrt(vals[1]))
        IF n_elements(where(times GT cl OR times LT cll)) GT 1 THEN BEGIN
          
        ENDIF ELSE BEGIN 
           
               dft=fft(times,-1) 
               dft=temporary(dft[0:nt/2-1])
               phold=[[phold],[2.*abs(dft)^2]]

        ENDELSE

     ENDIF 
     
  ENDFOR
  rpow_elem=[[rpow_elem],[n_elements(phold[0,*])]]
  pelm=n_elements(phold[0,*])
   
 
  FOR jj=0,nt/2-1 DO BEGIN
      med=median(phold[jj,1:pelm-1])
      sd=sqrt((moment(phold[jj,1:pelm-1]))[1])
      ravpow[jj,i-1,0]=med
      ravpow[jj,i-1,1]=sd
     
  ENDFOR
 
ENDFOR

;Phi averaging
phi_ang=45
no_phi=360/phi_ang
phiavpow=fltarr(nt/2,no_phi)
phipow_elem=0.


IF keyword_set(do_phi) THEN BEGIN
FOR i=1,no_phi DO BEGIN

     
     wedge=where(phi gt phi_ang*(i-1) and phi lt phi_ang*i and r lt upper_r)
 ;print,phi_ang*(i-1),phi_ang*(i)
 ;print,n_elements(wedge)
 
     xc=wedge mod nx
     yc=wedge/nx
     
     phold=fltarr(nt/2)
     FOR j=0,n_elements(xc)-1 DO BEGIN
  
     no_zero=n_elements(where(apocube[xc(j),yc(j),*] eq 0.))
      
     IF no_zero lt 5 THEN BEGIN 
        dft=fft(apocube[xc(j),yc(j),*],-1)  
        dft=temporary(dft[0:nt/2-1]) 
       
        phold=[[phold],[2.*abs(dft)^2]]
        ;IF i EQ 50 THEN BEGIN
        ;   plot,apocube[xc(j),yc(j),*]
        ;   pause
        ;ENDIF
     ENDIF 
     
  ENDFOR
  ;print,n_elements(phold[0,*])
  phiavpow[*,i-1]=total(phold[*,1:n_elements(phold[0,*])-1],2)/n_elements(phold[0,*])


ENDFOR
ENDIF

END