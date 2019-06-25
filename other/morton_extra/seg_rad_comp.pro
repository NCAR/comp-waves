;PURPOSE: - Calculate radially averaged power spectra for wedge of CoMP data
;
;INPUTS:- index - initial of observations (only used to calculate the radial coordinates) 
;         date  - date of aligned cube_ivw
;
;OPTIONAL INPUTS: - maxoffset - maximum offset from cross-correlation in wave_tracking.pro
;                 - apod - degree of apodisation (<1.) - default is 0.1
;
;                   extra - string relating to additional characters in cube_ivw name
;                   dopwid - calculates average using Doppler width
;                   inten - calculates average using intensity 
;
;OUTPUTS: - ravpow - radially averaged power array of frequency vs height
;           phipow - angular-ly averaged power array of frequency vs angle
;           rplot  - array of radial coordinates in arcsec
;
;
;RESTRICTIONS: ANGLE SELECTION - only looks for values phi < angles <phi2
;                                where phi<phi2, needs to be able to isolate values
;                                between phi2 > and  <phi, e.g between 370 and 10 
;                                -IS NOT FIXED



;Cacluate r,phi values for image using coord_cart_heilo
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

  ; temporal detrending per px
  ;ttrend = fltarr(nt)
  ;tf = findgen(nt) + 1
  FOR ix=0, nx-1 DO BEGIN
    FOR iy=0, ny-1 DO BEGIN
        ts = cube[ix,iy,0:nt-1]
        ;res=poly_fit(findgen(nt),ts,1,yfit=fit)
        fit=mean(ts)
        apocube[ix,iy,*]=(ts-fit)*apodt
        
    ENDFOR
  ENDFOR
  
  ; done
  return,apocube
END






;----------------------------------------------------------------------------
PRO seg_rad_comp,date,index=index,ravpow=ravpow,phi=phi,maxoffset=maxoffset,$           
                            rplot=rplot,apod=apod,extra=extra,dopwid=dopwid,inten=inten,$
                            outcube=outcube,man_ang=man_ang,errors=errors,for_log=for_log


; paths where the CoMP data is located and
; where the .sav-files should be written to
inpath  = 'data/CoMP/wave_tracking_output/'+date+'/'
outpath = 'data/CoMP/wave_tracking_output/'+date+'/'

inpath_ind  = 'data/CoMP/'+strmid(date,8,4,/reverse)+'/'+strmid(date,3,2,/reverse)+'/'+strmid(date,1,2,/reverse)+'/'
files=find_files('*.fts.gz',inpath_ind)
mreadfits,files(0),index

IF NOT keyword_set(extra) THEN restore,inpath+'cube_ivw_'+date+'.sav' $
ELSE restore,inpath+'cube_ivw_'+date+'_'+extra+'.sav'

restore,inpath+'wave_angle_'+date+'_3_3.5mHz.sav' 

IF keyword_set(errors) THEN restore,inpath+'cube_ivw_errors_incsee.sav'


data=cube_v
IF keyword_set(dopwid) THEN data=cube_w
IF keyword_set(inten) THEN data=cube_i



sz=size(data)
nx=sz(1)
ny=sz(2)
nt=sz(3)
xscale=index.cdelt1 ;pixel size in arcsec


IF keyword_set(errors) THEN data_err=dv ELSE data_err=fltarr(nx,ny,nt)


; Upper and lower limit for the radius in
; !be aware that the maximum offset of the cross-correlation
; (if performed) is added to this value!
lower_r = 228.*xscale
upper_r = 288.*xscale

IF n_elements(maxoffset) LT 1 THEN maxoffset=0.
IF n_elements(apod) LT 1 THEN apod=0.1


image=data[*,*,0]

r=frad_coord(index,image,phi=phi)

; create mask
; mask out pixels to invert
mask=intarr(nx,ny)

good=where(r ge lower_r+maxoffset and r le upper_r)
mask(good)=1.

data_mask=data
FOR i=0,nt-1 DO data_mask[*,*,i]=data[*,*,i]*mask

;SELECT WEDGE TO AVERAGE OVER
set_plot,'x'
aia_lct, wave=193, /load
tvim, reform(cube_i[*,*,0]) > 0 < 25

window,1
loadct,4,/silent
tvim,wave_angle

IF NOT keyword_set(man_ang) THEN BEGIN

   message,/cont,'You will be asked to pick 2 points - work clockwise!!'
   message,/cont,'Pick wedge side 1'
   cursor,x1,y1,/down
   x1=fix(x1)
   y1=fix(y1)
   print,phi[x1,y1]
   phione=phi[x1,y1]

   message,/cont,'Pick wedge side 2'
   cursor,x2,y2,/down
   x2=fix(x2)
   y2=fix(y2)
   print,phi[x2,y2]
   phitwo=phi[x2,y2]
ENDIF ELSE BEGIN
   answs=1d
   read, answs, prompt='ENTER ANGLE PHI 1: '  
   phione=answs

   read, answs, prompt='ENTER ANGLE PHI 2: '  
   phitwo=answs

ENDELSE

contour,phi,levels=phione,/over
contour,phi,levels=phitwo,/over
plots,[nx/2,nx/2],[ny/2,ny],color=0

;IF phi2

no_r=floor((upper_r-lower_r-maxoffset)/xscale)
rplot=findgen(no_r)*xscale+lower_r+maxoffset

ravpow=fltarr(nt/2,no_r,3)
rpow_elem=0.

print,'Starting'

;Apodise the time-series for each pixel
apocube=apod3dcube(data_mask,apod)
apocube_err=apod3dcube(data_err,apod)

outcube=fltarr(nt)

;Radial averaging
FOR i=1,no_r DO BEGIN

  IF phione LT phitwo THEN BEGIN
     ring=where(r GT (lower_r+maxoffset+(i-1)*xscale) $
     AND r LT (lower_r+maxoffset+i*xscale) AND phi GT phione AND phi LT phitwo )
  ENDIF ELSE BEGIN
       ring=where(r GT (lower_r+maxoffset+(i-1)*xscale) AND r LT (lower_r+maxoffset+i*xscale) AND (phi GT phione))
       ring2=where(r GT (lower_r+maxoffset+(i-1)*xscale) AND r LT (lower_r+maxoffset+i*xscale) AND (phi LT phitwo) )
       
       ring=[ring,ring2]
         
  ENDELSE

 ;print,lower_r+maxoffset+(i-1)*xscale,lower_r+maxoffset+(i)*xscale
 ;print,n_elements(ring)
  xc=ring mod nx
  yc=ring/nx

  phold=fltarr((nt-1)/2+1)
  phold_err=fltarr((nt-1)/2+1)
  FOR j=0,n_elements(xc)-1 DO BEGIN
  
     no_zero=n_elements(where(apocube[xc(j),yc(j),*] eq 0.))
      
     IF no_zero lt 5 THEN BEGIN 
         times=reform(apocube[xc(j),yc(j),*])
         IF keyword_set(errors) THEN timeserr=reform(apocube_err[xc(j),yc(j),*])
   
        outcube=[[outcube],[reform(times)]]
        ;Excludes pixels with time-series that have spikes greater
        ;than pm 4 sigma from the mean
        
   
        vals=moment(times)
        cl=(vals[0]+4.*sqrt(vals[1]))
        cll=(vals[0]-4.*sqrt(vals[1]))
        IF n_elements(where(times GT cl OR times LT cll)) GT 1 THEN BEGIN
          
        ENDIF ELSE BEGIN 
           
               dft=fft(times,-1) 
               dft=temporary(dft[0:(nt-1)/2])
               phold=[[phold],[2.*abs(dft)^2]]

               IF keyword_set(errors) THEN BEGIN
                  dft=fft(timeserr,-1) 
                  dft=temporary(dft[0:(nt-1)/2])
                  phold_err=[[phold_err],[2.*abs(dft)^2]]
               ENDIF

        ENDELSE
        
     ENDIF 
     
  ENDFOR
  
  rpow_elem=[[rpow_elem],[n_elements(phold[0,*])]]
  pelm=n_elements(phold[0,*]) 


  IF pelm gt 2 THEN BEGIN
    
    
  FOR jj=0,(nt-1)/2-1 DO BEGIN

      
      IF keyword_set(for_log) THEN BEGIN
            ;This setting calculates the weighted mean to log values
            ;ideal for log-log fits to data points     
            meanerr,alog(phold[jj,1:pelm-1]),phold_err[jj,1:pelm-1]/phold[jj,1:pelm-1],xmean,sigmam,sigmad
            ravpow[jj,i-1,0]=xmean
            ravpow[jj,i-1,1]=sigmad
            ravpow[jj,i-1,2]=pelm-2
            
      ENDIF ELSE BEGIN
            med=median(phold[jj,1:pelm-1])
            sd=sqrt((moment(phold[jj,1:pelm-1]))[1])
            ravpow[jj,i-1,0]=med
            ravpow[jj,i-1,1]=sd
            ravpow[jj,i-1,2]=pelm-2
      ENDELSE
  ENDFOR
  ENDIF

ENDFOR



END