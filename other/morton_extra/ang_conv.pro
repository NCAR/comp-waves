;convert CoMP angles to angles from normal to surface

PRO ang_conv, angles,index,oangles

oangles=angles
; Create arrays of r & phi values
r=frad_coord(index,angles,phi=phi)

in=where(phi lt 180)
oangles[in]=abs(phi[in]-90)+signum(phi[in]-90)*angles[in]

in=where(phi gt 180)
beta_a=(360-phi)
oangles[in]=abs(beta_a[in]-90)-signum(beta_a[in]-90)*angles[in]

oangles[where(index.mask eq 0)]=0
oangles[where(oangles gt 90)]-=180


END