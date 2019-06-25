
; duty_lim - lower limit on duty cycle
; frame_lim - lower limit on frame number
;

pro good_waves,duty_lim,frame_lim

fname='analysis/comp/good_waves.txt'
fmt='A,L,L,A,L,L,F,A,L,A,A,F' 
readcol,fname,start,date,time,ende,datee,timee,hour,hours,gi,good,dtuty,dc,format=fmt,/quick,count=count

in=where(dc gt duty_lim)
in2=where(gi gt frame_lim)
intersec=setintersection(in,in2)
print,'      Date     ',' Good frames ',' Duty cycle '
for i=0,n_elements(intersec)-1 do print,date[intersec[i]],gi[intersec[i]],dc[intersec[i]]

END