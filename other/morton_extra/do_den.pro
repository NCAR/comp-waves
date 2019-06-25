PRO do_den,den_map,lr_map

files=find_files('*.fts','data/comp/2012/03/27/new_flat_data/')
mreadfits,files(0),index
image=fltarr(620,620,3)
image2=fltarr(620,620,3)

FOR i=0,2 do image[*,*,i]=readfits(files(0),ext=i+1)
FOR i=0,2 do image2[*,*,i]=readfits(files(1),ext=i+1)

analytic_gauss_fit_fast,image[*,*,0],image[*,*,1],image[*,*,2],33.0,cube_v,cube_w,cube_i,/get_rid
analytic_gauss_fit_fast,image2[*,*,0],image2[*,*,1],image2[*,*,2],33.0,cube_v2,cube_w2,cube_i2,/get_rid

sz=size(image)
nx=sz(1)
ny=sz(2)
xscale=index(0).cdelt1 ;pixel size in arcsec
lower_r = 228.*xscale
upper_r = 288.*xscale
IF n_elements(maxoffset) LT 1 THEN maxoffset=0. ELSE maxoffset=maxoffset*xscale

;Calculate maps or r and phi values
r=frad_coord(index(0),image,phi=phi)

; create mask
; mask out pixels to invert
mask=intarr(nx,ny)

good=where(r ge lower_r+maxoffset and r le upper_r)
mask(good)=1.

lr_map=fltarr(nx,ny)
den_map=fltarr(nx,ny)

no_r=floor((upper_r-lower_r-maxoffset)/xscale)
rplot=findgen(no_r)*xscale+lower_r+maxoffset

restore,'data/comp/line_ratio/rat.sav'


FOR i=1,no_r-1 DO BEGIN                                                                      
    ring=where(r gt (lower_r+maxoffset+(i-1)*xscale) and r lt (lower_r+maxoffset+i*xscale) )
    xc=ring mod nx                                                                          
    yc=ring/nx                                                                              
   FOR j=0,n_elements(xc)-1 DO BEGIN                                                        
      no_zero=(where(cube_i[xc(j),yc(j)] eq 0.))                                      
      no_zero2=(where(cube_i2[xc(j),yc(j)] eq 0.))                                    
  IF no_zero[0] EQ -1 THEN no_zero=0 ELSE no_zero=n_elements(no_zero)                       
      IF no_zero2[0] EQ -1 THEN no_zero2=0 ELSE no_zero2=n_elements(no_zero2)               
      IF (no_zero lt 1) AND (no_zero2 lt 1) THEN BEGIN 
 av_int1074=(reform(cube_i[xc(j),yc(j)]))  
 av_int1079=(reform(cube_i2[xc(j),yc(j)]))                                            
 ratio=av_int1079/av_int1074 
         dense=calc_den(ratio_array,h,den,ratio,rplot(i))
         lr_map[xc(j),yc(j)]=ratio
         den_map[xc(j),yc(j)]=dense  
            
      ENDIF 
      
   ENDFOR
ENDFOR
 
end

