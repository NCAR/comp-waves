pro test_memcof
 
;  generate test function with noise
 
;restore,'data/comp/wave_tracking_output/20120411/cube_ivw_20120411.sav',/verb
;data=reform(cube_v[152,124,50:150])
restore,'analysis/comp/mem_test.sav'
data=data[0:199]
data=data-mean(data) 
;data=smooth(data,3,/edge_truncate)

;freq1= (3.3 + randomn(seed,1))/16.666
;freq2= (3.0 + randomn(seed,1))/16.666
;freq3= (2.7 + randomn(seed,1))/16.666
;data=0.33333*sin(findgen(200)*2.*!pi*freq1(0)) + $
;     0.33333*sin(findgen(200)*2.*!pi*freq2(0)) + $
;     0.33333*sin(findgen(200)*2.*!pi*freq3(0)) + $
;     randomn(seed,200)/20.
orig_data=data
 
;  create gaps
 
gap_start=intarr(3)
gap_end=intarr(3)
 
gap_start(0)=6
gap_end(0)=10
gap_start(1)=105
gap_end(1)=109
gap_start(2)=195
gap_end(2)=198
ngaps=3
 
for i=0,ngaps-1 do data(gap_start(i):gap_end(i))=0.
 
;  plot data
 
plot, orig_data, xrange=[0,200], xstyle=1


;b=memcof(data[6:40], 35, xms, coef)
;fixrts, coef, 35
;future=predic(data[6:40],coef,40)
 

;stop 
;  fill gaps
fill_gaps, data, gap_start, gap_end, ngaps,coef_num=30,/verb,/aic_sel
 
for i=0,ngaps-1 do begin
  gap_size=gap_end(i)-gap_start(i)+1
  oplot, findgen(gap_size)+gap_start(i), data(gap_start(i):gap_end(i)), psym=1
endfor

stop
 
;ans=' '
;read,'enter return',ans
diff=orig_data-data
plot, diff, xrange=[0,200], xstyle=1

end
