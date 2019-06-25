;
;
;
;INPUTS: x1,x2 - heights to fit between
;

PRO fit_slope,data,x1=x1,x2=x2,xscale=xscale,fitres=fitres

sz=size(data)
no_freq=sz(1)
height=sz(2)

IF n_elements(x1) LT 1 THEN x1=0
IF n_elements(x2) LT 1 THEN x2=height-1
IF n_elements(xscale) LT 1 THEN xscale=4.46*0.725 ; CoMP pix size in Mm
IF n_elements(dt) LT 1 THEN dt=30. ; standard CoMP cadence
f=findgen(no_freq)/(2.*no_freq*dt)

x=findgen(x2-x1+1)*xscale

;Provide an initial guess of the function’s parameters.


;Can we weight the function?
weights=fltarr(x2-x1+1)
weights[*]=1.

fitres=fltarr(no_freq,4)

FOR i=0,no_freq-1 DO BEGIN
    
    plot,data[i,x1:x2]
    dat=reform(data[i,x1:x2])
 
    ;Compute the parameters.
    start=[max(dat),-0.1]
    result = mpfitfun('myexp',X,dat,weights,start,perror=sigma,/quiet)

    oplot,result[0]*exp(x*result[1])
    ;print,result
    fitres[i,0:1]=result
    fitres[i,2:3]=sigma
    pause
   
ENDFOR
plot_oo,f,abs(1./fitres[*,1]),xrange=[1e-4,0.1],yrange=[1,100]
res=poly_fit(alog10(f(1:15)),alog10(abs(1/fitres[1:15,1])),1)
print,res
oplot,f,10^(res[0])*f^(res[1])
END


