function parabola_test,x,y

;  function to return minimum position of parabola given 3 data points
f=x[1]-(y[2]-y[1])/(y[2]-2.*y[1]+y[0])-0.5

res=poly_fit(x,y,2,sigma=sigma,yfit=yfit)
 print,-res[1]/res[2]/2.
return,f

end
