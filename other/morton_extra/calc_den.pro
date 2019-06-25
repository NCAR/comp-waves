;Calculate density from 2D grid of line ratio vs height
;
FUNCTION calc_den,rata,h,den,ratio,height,index=index

  ;height input in arcsec, convert to solar radii
  ;xscale in Mm
  rsun_mm=695.5
  height=height/index.cdelt1*index.xscale/rsun_mm

  hloc=where(h lt height+0.01 and h gt height-0.01)
  
  IF n_elements(hloc) LT 3 THEN hloc=[hloc(0)-1,hloc(0),hloc(0)+1] 
  hn=n_elements(hloc)
  
  ;The pm0.05 range should cope with the steepest part of the ratio slope
  ; and provide two values to interpolate between
  in=rata[hloc(1),*]
  
  loc=where((in GT (ratio[0]-0.05)) AND (in LT (ratio[0]+0.05)))
  ;print,loc,ratio+0.05,ratio-0.05

  IF n_elements(loc) LT 3 THEN BEGIN
     IF n_elements(loc) LT 2 THEN loc=[loc(0)-1,loc(0),loc(0)+1]$
     ELSE loc=[loc(0),loc(1),loc(1)+1] 
  ENDIF
  ln=n_elements(loc)
  
 
  IF (loc[0] GT -1) AND (loc[2] lt n_elements(rata[0,*])) THEN BEGIN
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
  ENDIF ELSE dense=0

return,dense

END