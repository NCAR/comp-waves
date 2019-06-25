function memcof,data,m,xms,d_j

;  memcof subroutine from Numerical Methods, 
;  Based on Burg, J.P. 1968, “A New Analysis Technique for Time Series Data,” reprinted in Childers, 1978.
;
;  INPUTS - data
;         - m - number of Linear Prediction (LP) coefficients required
;
;  OUTPUTS - xms - mean square discrepancy
;          - d_j - Estimated LP coefficients


;  given a real vector of DATA of length n, and given M, this routine
;  returns a vector d_j of length M with d_j(j)=a_j, and a scalar PM=b_0,
;  which are the coefficients for Maximum Entropy Method spectral estimation.
;  The user must provide workspace vectors WK1, WK2, and WKM of lengths
;  N, N, and M, respectively.

n=n_elements(data)
nmax=10000 & mmax=200
if m gt mmax or n gt nmax then return,-1

d_j=fltarr(m)
wk1=fltarr(n)
wk2=fltarr(n)
wkm=fltarr(m)

xms=total(data^2)/float(n)

wk1[0]=data[0:n-2]
wk2[0]=data[1:n-1]


for k=1,m do begin
  pneum=total(wk1[0:n-k-1]*wk2[0:n-k-1])
  denom=total(wk1[0:n-k-1]^2+wk2[0:n-k-1]^2)
  
  d_j[k-1]=2.*pneum/denom
  xms=xms*(1.-d_j[k-1]^2)

  for i=1,k-1 do d_j[i-1]=wkm[i-1]-d_j[k-1]*wkm[k-i-1]

  if k eq m then return,0

  wkm[0:k-1]=d_j[0:k-1]

  for j=1,n-k-1 do begin
    wk1[j-1]=wk1[j-1]-wkm[k-1]*wk2[j-1]
    wk2[j-1]=wk2[j]-wkm[k-1]*wk1[j]
  endfor
endfor

return,-1
end
