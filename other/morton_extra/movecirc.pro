

    cube_v=fltarr(620,620,1024)
    
    veloc=fltarr(270,320,1024)
    veloc[20:269,0:319,*]=velocity
    for i=0, 1023 do begin
    
    for j=0,249 do begin
    	for k=0,319 do begin
    		cube_v[j+9,k+41,i]=veloc[j,k,i]
    	endfor
    endfor
    
    
    endfor

nx=620
ny=620
nt=1024
xscale=0.6
yscale=0.6
norm_cadence=29.0



index=create_struct( 'NAXIS1', nx, 'NAXIS2', ny, 'NFRAMES', nt, $ 
                    'XSCALE', xscale, 'NORM_CADENCE',  norm_cadence)

print,'pick points on occulting disk edge'
window,xs=1000,ys=1000
plot_image,cube_v(*,*,0)
pick,x1,y1
pick,x2,y2
pick,x3,y3
xn=[x1,x2,x3]
yn=[y1,y2,y3]
cir_3pnt,xn,yn,r0,x0,y0
print,'radius is ',r0 
print,'center is', x0,y0
;;;;;;;;;;;;;;;;;;;;;;

;;;;;; defining mask for velocity images;;;;;
lower_r=r0+5
upper_r=lower_r+80.




print,'radius is ',r0 
print,'center is', x0,y0
;;;;;;;;;;;;;;;;;;;;;;

;;;;;; defining mask for velocity images;;;;;
lower_r=r0+5
upper_r=lower_r+80.
save,r0,x0,y0,filename='circle.sav'

index=create_struct( 'NAXIS1', nx, 'NAXIS2', ny, 'NFRAMES', nt, $ 
                    'XSCALE', xscale, 'NORM_CADENCE',  norm_cadence,'LOWER_R',lower_r,'UPPER_R',upper_r)

save,index,cube_v,filename='cube_20051030.sav'
end