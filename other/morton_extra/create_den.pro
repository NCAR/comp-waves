;PURPOSE - Calculate the density from ratio of CoMP lines.
;          Calculates the ratio of the 10747 and 10798 lines and
;          uses a pre-computed look-up table of intensity ratios
;          created with CHIANTI for a set T and photo-ionisation rate.
;

PRO create_den,inpath,den_map,lr_map

restore,inpath+'aligned7479.sav'

sz=size(cubei1074)
nx=sz(1)
ny=sz(2)
IF sz(0) EQ 3 THEN nt=sz(3) ELSE nt=1

xscale=index.cdelt1 ;pixel size in arcsec
lower_r = 225*xscale
upper_r = 280*xscale
maxoffset=0

;Calculate maps or r and phi values
;frad_coord calculates values using date_obs in UTC
;CoMP index.date_obs default is in local time - Manua Loa
;However, format causes problems with coord_cart_helio
index.date_obs=index.date_d$obs+' '+index.time_d$obs
r=frad_coord(index,cubei1074[*,*,0],phi=phi)

; create mask
; mask out pixels to invert
;mask=index.mask
mask=fltarr(nx,ny)
mask[where(r gt lower_r and r lt upper_r)]=1.

lr_map=fltarr(nx,ny,nt)
den_map=fltarr(nx,ny,nt)

no_r=floor((upper_r-lower_r-maxoffset)/xscale)
rplot=findgen(no_r)*xscale+lower_r+maxoffset

restore,'analysis/comp/line_ratio/rat.sav'


;For different values of r
FOR i=1,no_r-1 DO BEGIN                                                                      
    ring=where(r gt (lower_r+maxoffset+(i-1)*xscale) and r lt (lower_r+maxoffset+i*xscale) )
    xc=ring mod nx                                                                          
    yc=ring/nx

   ;For the different time-frames
   FOR kk=0,nt-1 DO BEGIN

          cube_i=reform(cubei1074[*,*,kk])
          cube_i2=reform(cubei1079[*,*,kk])                                                                     
   
           FOR j=0,n_elements(xc)-1 DO BEGIN                                                        
                
                no_zero=(where(cube_i[xc(j),yc(j)] eq 0.))                                      
                no_zero2=(where(cube_i2[xc(j),yc(j)] eq 0.))                                    
                IF no_zero[0] EQ -1 THEN no_zero=0 ELSE no_zero=n_elements(no_zero)                       
                IF no_zero2[0] EQ -1 THEN no_zero2=0 ELSE no_zero2=n_elements(no_zero2)               
               
                IF (no_zero lt 1) AND (no_zero2 lt 1) THEN BEGIN 
                       av_int1074=(reform(cube_i[xc(j),yc(j)]))  
                       av_int1079=(reform(cube_i2[xc(j),yc(j)]))                                            
                       ratio=av_int1079/av_int1074 
                       dense=calc_den(ratio_array,h,den,ratio,rplot(i),index=index)
                       lr_map[xc(j),yc(j),kk]=ratio
                       den_map[xc(j),yc(j),kk]=dense  
                ENDIF 
           ENDFOR
   ENDFOR
ENDFOR
 
end

