;PURPOSE - test correlations between i, v and w
;
;
;
;
;



 ; apodizes time-series, with optional de-trending 
;----------------------------------------------------------------------------
FUNCTION apod_ts,ts,apod ;,detrend=detrend
    
  ; get cube dimensions
  sizecube=size(ts)
  nt=sizecube[1]
  apocube = fltarr(nt)

  ; define temporal apodization
  apodt = fltarr(nt)+1
  IF (apod ne 0) THEN BEGIN
    apodrimt = nt*apod
    apodt[0] = (sin(!pi/2.*findgen(apodrimt)/apodrimt))^2
    apodt = apodt*shift(rotate(apodt,2),1)   
  ENDIF

  ; temporal detrending, not per px, only mean-image trend 
  
  res=poly_fit(findgen(nt),ts,1,yfit=ttrend)

  outim=ts
  outim= (ts-ttrend)*apodt
  apocube=outim 

  
  ; done
  return,apocube
END




PRO cc_ivw,cc=cc

restore,'data/comp/wave_tracking_output/20120327/cube_ivw_20120327.sav'

apod=0.1


sz=size(cube_i)
nx=sz(1)
ny=sz(2)
nt=sz(3)

cc=fltarr(nx,ny,6)


FOR i=0,nx-1 DO BEGIN
     FOR j=0,ny-1 DO BEGIN
    
     tsv=reform(cube_v[i,j,*])
     tsw=reform(cube_w[i,j,*])
     tsi=reform(cube_i[i,j,*])

     tsva=apod_ts(tsv,apod)
     tswa=apod_ts(tsw,apod)
     tsia=apod_ts(tsi,apod)

     cc[i,j,0]=max(c_correlate(tsva,tswa,findgen(nt)-nt/2),loc)
     cc[i,j,3]=loc
     cc[i,j,1]=max(c_correlate(tsva,tsia,findgen(nt)-nt/2),loc)
     cc[i,j,4]=loc
     cc[i,j,2]=max(c_correlate(tsia,tswa,findgen(nt)-nt/2),loc)
     cc[i,j,5]=loc

  
     counter,i*ny+j,nx*ny,/percent
     ENDFOR
     

ENDFOR

 im=cc[*,*,0]
 im2=cc[*,*,3]
 im[where(finite(im,/nan,sign=-1))]=0D
 im[where(finite(im,/nan,sign=-1))]=0D
 im2[where(finite(im,/nan,sign=1))]=0D
 im2[where(finite(im,/nan,sign=1))]=0D
 cc[*,*,0]=im
 cc[*,*,3]=im2

 im=cc[*,*,1]
 im2=cc[*,*,4]
 im[where(finite(im,/nan,sign=-1))]=0D
 im[where(finite(im,/nan,sign=1))]=0D
 im2[where(finite(im,/nan,sign=1))]=0D
 im2[where(finite(im,/nan,sign=1))]=0D
 cc[*,*,1]=im
 cc[*,*,4]=im2

 im=cc[*,*,2]
 im2=cc[*,*,5]
 im[where(finite(im,/nan,sign=-1))]=0D
 im[where(finite(im,/nan,sign=1))]=0D
 im2[where(finite(im,/nan,sign=1))]=0D
 im2[where(finite(im,/nan,sign=1))]=0D
 cc[*,*,2]=im
 cc[*,*,5]=im2


END