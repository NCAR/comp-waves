;PURPOSE: - Calculate radially and angular-ly averaged power spectra for CoMP data
;
;INPUTS:- index - initial of observations (only used to calculate the radial coordinates) 
;         date  - date of aligned cube_ivw
;
;OPTIONAL INPUTS: - maxoffset - maximum offset from cross-correlation in wave_tracking.pro
;                 - apod - degree of apodisation (<1.) - default is 0.1
;
;                   extra - string relating to additional characters in cube_ivw name
;                   
;                   /do_phi - will provide the data averaged over different phi angles (no radial seperation)
;
;
;OUTPUTS: - av_den - radially averaged density vs height
;           phipow - angular-ly averaged power array of frequency vs angle
;           rplot  - array of radial coordinates in arcsec
;
;



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

FUNCTION calc_den,rata,h,den,ratio,height

  ;height input in arcsec, convert to solar radii
  height=height*725./6.955e5
  hloc=where(h lt height+0.01 and h gt height-0.01)
  
  IF n_elements(hloc) LT 3 THEN hloc=[hloc(0)-1,hloc(0),hloc(0)+1] 
  hn=n_elements(hloc)
  
  ;The pm0.05 range should cope with the steepest part of the ratio slope
  ; and provide two values to interpolate between
  loc=where((rata[hloc(1),*] LT (ratio+0.05)) AND (rata[hloc(1),*] GT (ratio-0.05)))
 ; plot,rata[pos,*]
 ; plots,[0,80],[ratio+0.05,ratio+0.05]
 ; plots,[0,80],[ratio-0.05,ratio-0.05]
 

  IF n_elements(loc) LT 3 THEN BEGIN
     IF n_elements(loc) LT 2 THEN loc=[loc(0)-1,loc(0),loc(0)+1]$
     ELSE loc=[loc(0),loc(1),loc(1)+1] 
  ENDIF
  ln=n_elements(loc)

  
  ;2D fit to surface defined by height(x) vs density(y) and ratio (z), i.e. z=f(x,y)
  res=sfit(rata[hloc(0):hloc(hn-1),loc(0):loc(ln-1)],2,kx=kx)    

  ;ratio=f(x)*y^2+g(x)*y+h(x)
  ;We know ratio and height need to find y
  ;Use quadratic formulae to solve for y and choose gradient gt 0
  
  x=[1,height,height^2]
  a=total(kx[0,*]*x)             
  b=total(kx[1,*]*x)             
  c=total(kx[2,*]*x)-ratio       
  m1=(-b+sqrt(b^2-4*a*c))/2./a
  m2=(-b-sqrt(b^2-4*a*c))/2./a
  m=m1>m2
   
  den_scale=(den-shift(den,1))[1]  
  dense=den(loc(0))+m*den_scale

return,dense

END


;----------------------------------------------------------------------------
PRO rad_ave_comp_den,date,index=index,av_den=av_den,phiavpow=phiavpow,maxoffset=maxoffset,$           
                            rplot=rplot,nfiles=nfiles,do_phi=do_phi


; paths where the CoMP data is located and
; where the .sav-files should be written to

outpath = 'data/CoMP/wave_tracking_output/'+date+'/'

inpath  = 'data/CoMP/'+strmid(date,8,4,/reverse)+'/'+strmid(date,3,2,/reverse)+'/'+strmid(date,1,2,/reverse)+'/'+date+'.comp.1074.daily_dynamics.5/'
files=find_files('*.fts.gz',inpath)

test=readfits(files(0),ext=1,/silent)
sz=size(test)
data1074=fltarr(sz(1),sz(2),nfiles)
data1074[*,*,0]=test
FOR i=1,nfiles-1 DO data1074[*,*,i]=readfits(files(i),ext=1,/silent)

inpath  = 'data/CoMP/'+strmid(date,8,4,/reverse)+'/'+strmid(date,3,2,/reverse)+'/'+strmid(date,1,2,/reverse)+'/'+date+'.comp.1079.daily_dynamics.5/'
files=find_files('*.fts.gz',inpath)

test=readfits(files(0),ext=1,/silent)
sz=size(test)
data1079=fltarr(sz(1),sz(2),nfiles)
data1079[*,*,0]=test
FOR i=1,nfiles-1 DO data1079[*,*,i]=readfits(files(i),ext=1,/silent)

mreadfits,files(0),index


nx=sz(1)
ny=sz(2)
xscale=index(0).cdelt1 ;pixel size in arcsec


; Upper and lower limit for the radius in
; !be aware that the maximum offset of the cross-correlation
; (if performed) is added to this value!
lower_r = 228.*xscale
upper_r = 288.*xscale

IF n_elements(maxoffset) LT 1 THEN maxoffset=0. ELSE maxoffset=maxoffset*xscale



image=data1079[*,*,0]

r=frad_coord(index(0),image,phi=phi)

; create mask
; mask out pixels to invert
mask=intarr(nx,ny)

good=where(r ge lower_r+maxoffset and r le upper_r)
mask(good)=1.

data_mask1074=data1074
data_mask1079=data1079
for i=0,nfiles-1 do data_mask1074[*,*,i]=data1074[*,*,i]*mask
for i=0,nfiles-1 do data_mask1079[*,*,i]=data1079[*,*,i]*mask



no_r=floor((upper_r-lower_r-maxoffset)/xscale)
rplot=findgen(no_r)*xscale+lower_r+maxoffset

av_den=fltarr(no_r,3)

restore,'data/comp/line_ratio/rat.sav'

;Radial averaging
FOR i=1,no_r DO BEGIN

   ring=where(r gt (lower_r+maxoffset+(i-1)*xscale) and r lt (lower_r+maxoffset+i*xscale) )
   xc=ring mod nx
   yc=ring/nx

  phold=0.
  pholde=0.
  FOR j=0,n_elements(xc)-1 DO BEGIN
  
     no_zero=(where(data_mask1074[xc(j),yc(j),*] eq 0.))
     no_zero2=(where(data_mask1079[xc(j),yc(j),*] eq 0.))
     ;help,no_zero,no_zero2

     IF no_zero[0] EQ -1 THEN no_zero=0 ELSE no_zero=n_elements(no_zero)
     IF no_zero2[0] EQ -1 THEN no_zero2=0 ELSE no_zero2=n_elements(no_zero2)
      
     IF (no_zero lt 1) AND (no_zero2 lt 1) THEN BEGIN 
        av_int1074=total(reform(data_mask1074[xc(j),yc(j),*]))/nfiles
        av_int1079=total(reform(data_mask1079[xc(j),yc(j),*]))/nfiles
        
        ratio=av_int1079/av_int1074
        ;print,av_int1079,av_int1074,ratio,max(ratio_array)
        dense=calc_den(ratio_array,h,den,ratio,rplot(i))
        phold=[phold,dense]       
   
        ;Basic error calculation - variation of intensity over the images
        ;er_int1074=sqrt((moment(data_mask1074[xc(j),yc(j),*]))[1])
        ;er_int1079=sqrt((moment(data_mask1079[xc(j),yc(j),*]))[1]) 
        ;err_rat=(er_int1074/av_int1079)^2+(er_int1079*av_int1074/av_int1079^2)^2
        ;pholde=[pholde,sqrt(err_rat)]
        
     ENDIF 
     
  ENDFOR
  
  pelm=n_elements(phold[*])
  IF pelm gt 2 THEN BEGIN
     plothist,phold,bin=0.1
     av_den(i,0)=mean(phold[1:pelm-1])
     av_den(i,1)=sqrt((moment(phold[1:pelm-1]))[1])
     av_den(i,2)=pelm
     ;pause
  ENDIF
 
ENDFOR

;Phi averaging
;phi_ang=45
;no_phi=360/phi_ang
;phiavpow=fltarr(nt/2,no_phi)
;phipow_elem=0.


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