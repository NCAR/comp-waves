pro fixrts,d,npoles
 
;  subroutine fixrts from Numerical Recipes
 
;  Given the LP coefficients D(j), j=1,npoles, this routine finds all roots
;  of the characteristic polynomial, reflects any roots that are outside the
;  unit circle back inside, and returns a modified set of d(j)'s.
 
npmax=npoles+1	;maximum number of poles
a=complexarr(npmax)
roots=complexarr(npmax)
 
a[npoles]=complex(1.,0.)
 
for j=npoles,1,-1 do a[j-1]=complex(-d[npoles-j],0.)
 
;Find roots of equation 
roots=fz_roots(a,/double)

;Look for roots outside unit circle and reflect back in 
index=where(abs(roots) gt 1)
roots[index]=1./conj(roots[index])
;for j=1,npoles do if(abs(roots[j-1]) gt 1.) then $
;  roots[j-1]=1./conj(roots[j-1])
 

; Reconstruct polynomial coefficients 
a[0]=-roots[0]
a[1]=complex(1.,0.)
 
for j=2,npoles do begin
  a[j]=complex(1.,0.)
  for i=j,2,-1 do a[i-1]=a[i-2]-roots[j-1]*a[i-1]
  a[0]=-roots[j-1]*a[0]
endfor
 
d=-reverse(float(a[0:-2]))
end
