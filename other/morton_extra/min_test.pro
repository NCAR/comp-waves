function min_test,x,y,error=error

xb=0.;mean(x)
yb=0.;mean(x)
sxx=total((x-xb)^2)
syy=total((y-yb)^2)
sxy=total((x-xb)*(y-yb))

if sxy ne 0. then begin

 ;compute two possible choices for slope

  arg1=(-sxx+syy)/(2.*sxy)
  arg2=sqrt( (sxx-syy)^2 +4.*sxy^2 )/(2.*sxy)
  m1=arg1-arg2
  m2=arg1+arg2

;  find sum of distance squared for two slopes and choose slope which minimizes sum of distance squared

  distsq1=(m1^2*sxx  -2.*m1*sxy +syy)/(1.+m1^2)
  distsq2=(m2^2*sxx  -2.*m2*sxy +syy)/(1.+m2^2)

  if distsq1 le distsq2 then m=m1 else m=m2
  b4=(m)
  

b1=sxy/sxx
varb1=1/sxx^2*total((x-xb)^2*(y-b1*x-yb+b1*xb)^2)
b2=syy/sxy
varb2=1/sxy^2*total((y-yb)^2*(y-b2*x-yb+b2*xb)^2)

;b4=0.5*( (b2-1./b1) + sign(1,sxy)* sqrt( 4+(b2-1./b1)^2 ) )
cov=1./(b1*sxx^2)*total( (x-xb)* (y-yb)* (y-mean(y)-b1*(x-xb))* (y-yb-b2*(x-xb)) )
varb4=b4^2/(4.*b1^2+(b1*b2-1)^2)*(1./b1^2*varb1+2.*cov+b1^2*varb2)

b4=atan(b4)
error=varb4

endif else if sxx ge syy then b4=0. else b4=!pi/2.

return,b4

END