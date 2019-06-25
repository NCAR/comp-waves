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
;                                -NEEDS FIXING - Logic statements don't work



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

;#############################################################

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
PRO seg_rad_comp_den,date,index=index,av_den=av_den,maxoffset=maxoffset,$           
                            rplot=rplot,nfiles=nfiles

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

inpath  = 'data/CoMP/wave_tracking_output/'+date+'/'
restore,inpath+'wave_angle_'+date+'_3_3.5mHz.sav' 

nx=sz(1)
ny=sz(2)
xscale=index.cdelt1 ;pixel size in arcsec


; Upper and lower limit for the radius in
; !be aware that the maximum offset of the cross-correlation
; (if performed) is added to this value!
lower_r = 228.*xscale
upper_r = 288.*xscale

IF n_elements(maxoffset) LT 1 THEN maxoffset=0.
IF n_elements(apod) LT 1 THEN apod=0.1


image=data1079[*,*,0]

r=frad_coord(index,image,phi=phi)

; create mask
; mask out pixels to invert
mask=intarr(nx,ny)

good=where(r ge lower_r+maxoffset and r le upper_r)
mask(good)=1.

data_mask1074=data1074
data_mask1079=data1079
for i=0,nfiles-1 do data_mask1074[*,*,i]=data1074[*,*,i]*mask
for i=0,nfiles-1 do data_mask1079[*,*,i]=data1079[*,*,i]*mask

;SELECT WEDGE TO AVERAGE OVER
set_plot,'x'
aia_lct, wave=193, /load
tvim, reform(data_mask1074[*,*,0]) > 0 < 25

window,1
loadct,4,/silent
tvim,wave_angle

;message,/cont,'You will be asked to pick 2 points - work clockwise!!'
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


contour,phi,levels=phione,/over
contour,phi,levels=phitwo,/over
plots,[nx/2,nx/2],[ny/2,ny],color=0

;IF phi2

no_r=floor((upper_r-lower_r-maxoffset)/xscale)
rplot=findgen(no_r)*xscale+lower_r+maxoffset

av_den=fltarr(no_r,3)

print,'Starting'



restore,'data/comp/line_ratio/rat.sav'

;Radial averaging
FOR i=1,no_r DO BEGIN

  IF phione LT phitwo THEN BEGIN
     ring=where(r GT (lower_r+maxoffset+(i-1)*xscale) $
     AND r LT (lower_r+maxoffset+i*xscale) AND phi GT phione AND phi LT phitwo )
  ENDIF ELSE BEGIN
       ring=where(r GT (lower_r+maxoffset+(i-1)*xscale) AND r LT (lower_r+maxoffset+i*xscale) AND (phi GT phione))
       ring2=where(r GT (lower_r+maxoffset+(i-1)*xscale) AND r LT (lower_r+maxoffset+i*xscale) AND (phi LT phitwo) )
       ;help,ring,ring2
       ring=[ring,ring2]
         
  ENDELSE
 ;print,lower_r+maxoffset+(i-1)*xscale,lower_r+maxoffset+(i)*xscale
 ;print,n_elements(ring)
  xc=ring mod nx
  yc=ring/nx

  phold=0.
  pholde=0.
  print,n_elements(ring),i
  FOR j=0,n_elements(xc)-1 DO BEGIN
  
     no_zero=(where(data_mask1074[xc(j),yc(j),*] eq 0.))
     no_zero2=(where(data_mask1079[xc(j),yc(j),*] eq 0.))

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
        print,'exiin'
     ENDIF 
          
      
     
  ENDFOR
  
 
  pelm=n_elements(phold[0,*]) 


 pelm=n_elements(phold[*])
  IF pelm gt 2 THEN BEGIN
     plothist,phold,bin=0.1
     av_den(i,0)=mean(phold[1:pelm-1])
     av_den(i,1)=sqrt((moment(phold[1:pelm-1]))[1])
     av_den(i,2)=pelm
     ;pause
  ENDIF

ENDFOR



END