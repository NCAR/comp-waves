FUNCTION predic, data, d_j, nfut
 
;  INPUTS - data values
;         - LP coefficients
;         - number of coefficients to predict
;   
;  subroutine predic from Numerical Recipes
 
;  Given DATA(j), j=1,ndata, and given the data's LP coefficients D(i),
;  i=1,npoles, this routine predicts the next nfut data points, which
;  it returns in the array future. Note that the routine references only
;  the last npoles values of data, as initial values for the prediction.
 
future=fltarr(nfut)
;npmax=100	;maximum number of allowable poles
np=n_elements(d_j)
reg=fltarr(np)
 
reg=reverse(data)
 
for j=1,nfut do begin
  sum=total( d_j*reg )
  reg=shift(reg,1)
  reg(0)=sum    
  future(j-1)=sum
endfor

return,future 
end
