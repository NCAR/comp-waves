PRO loop_amplitude

;inpath='data/comp/wave_tracking_output/20120410/'
;outpath=inpath+'test/'
;loop=lp_read(inpath+'SPLINEcuts/velocity/rotTDR_intenscube_V1.fcube')
;restore,inpath+'20120410_loopoutput.sav'
;loop=fltarr(20,320,94)
;for i=0,93 do loop[*,*,i]=reform(output[*,i,1,*])

inpath='data/comp/wave_tracking_output/20130914/'
outpath=inpath
restore,inpath+'loop_profiles.sav'
loop=out_v[*,5:170,*]

vmap=transpose(loop,[1,2,0])

sz=size(vmap)
npt=sz(2) & nt=sz(1) & nslice=sz(3)


vmap=rebin(vmap,nt/2,npt,nslice)


stop

sz=size(vmap)
npt=sz(2) & nt=sz(1) & nslice=sz(3)
num_mc=500

freq=findgen(nt/2)/nt/30.

; define temporal apodization
apodt = fltarr(nt)+1
apod=0.1
apodrimt = nt*apod
apodt[0] = (sin(!pi/2.*findgen(apodrimt)/apodrimt))^2
apodt = apodt*shift(rotate(apodt,2),1) 
CPGT=total(apodt)/n_elements(apodt) ;Coherent Power Gain - correction factor needed for power after apodisation

FOR i=0,npt-1 DO FOR j=0,nslice-1 DO BEGIN
   
    res=poly_fit(findgen(nt),vmap(*,i,j),1,yfit=yfit)
    vmap(*,i,j)=vmap(*,i,j)-yfit ;remove temporal mean
   ;vmap(*,i,j)=vmap(*,i,j)-mean(vmap(*,i,j)) ;remove temporal mean
ENDFOR 


vh=fltarr(nt,npt,nslice,num_mc)
vh2=fltarr(nt,npt,nslice,num_mc)
vmap_pro=fltarr(nt,npt,nslice,num_mc)
vmap_ret=fltarr(nt,npt,nslice,num_mc)



rndser=randomn(systime,nt,npt,nslice,num_mc)
;estimate of noise on data
sds=fltarr(npt,nslice)

FOR j=0,npt-1 DO FOR k=0,nslice-1 DO sds[j,k]=sqrt(mean((reform(vmap[*,j,k])-smooth(reform(vmap[*,j,k]),3,/edge_truncate))^2)) 

FOR k=0, num_mc-1 DO BEGIN
FOR i=0,nslice-1 DO BEGIN

   err_arr=rebin(transpose(reform(sds[*,i])),nt,npt)
   dum=smooth(reform(vmap(*,*,i)),[3,1],/edge_truncate)+rndser[*,*,i,k]*err_arr
   
   trans=fft(gauss_smooth(dum,[0.5,0.5]),-1)       ;compute fourier transform

   ;mag=abs(trans)         ;compute magnitude and phase
   ;phase=atan(imaginary(trans),float(trans))
 
   ;back=rebin( rebin(mag(nt/2-(nt/2-2):nt/2+(nt/2-2),*),1,npt) ,nt,npt)
   ;mag=mag-back       ;remove high temporal frequency noise
 
   ;back=rebin(mag(*,npt/2),nt,npt)
   ;mag=mag-back       ;remove high spatial frequency noise
 
   ;mag=mag > 0.   ;insure that magnitude is positive
 
  ;pw_trans=complex(mag*cos(phase),mag*sin(phase))   ;recompute transform
   
  ;pw_trans=pw_trans*filter      ;remove low temporal frequencies;  filter is defined at top

  pro_trans=trans             ;select prograde waves (assume nt even, npt odd)
  pro_trans(1:nt/2-1,0:npt/2)=0.
  pro_trans(nt/2+1:nt-1,npt/2+1:npt-1)=0.
  pro_vel=float(fft(pro_trans,1))

  vmap_pro[*,*,i,k]=pro_vel

  ret_trans=trans             ;select retrograde waves
  ret_trans(1:nt/2-1,npt/2+1:npt-1)=0.
  ret_trans(nt/2+1:nt-1,0:npt/2)=0.
  ret_vel=float(fft(ret_trans,1))

  vmap_ret[*,*,i,k]=ret_vel

  FOR j=0,npt-1 DO vh[*,j,i,k]=2.*abs(fft(pro_vel[*,j]*apodt))/cpgt
  FOR j=0,npt-1 DO vh2[*,j,i,k]=2.*abs(fft(ret_vel[*,j]*apodt))/cpgt



ENDFOR
ENDFOR
up=(moment(mean(vh[0:nt/2,*,6:12,*],dim=3),dim=3))[*,*,0:1] 
down=(moment( mean(vh2[0:nt/2,*,6:12,*],dim=3) ,dim=3) )[*,*,0:1] 

TF=up[*,*,0]/rebin(reform(up[*,0,0]),nt/2+1,npt)


loadct,0
!p.background=255
window,0

plot,total(up[0:6,*,0],1)/(total(up[0:6,*,0],1))[3],color=0,/nodata
loadct,13  
oplot,smooth(total(up[1:3,*,0],1)/(total(up[1:3,*,0],1))[3],3),color=0             
oplot,smooth(total(up[4:6,*,0],1)/(total(up[4:6,*,0],1))[3],3)  ,linestyle=2,color=50  
oplot,smooth(total(up[7:9,*,0],1)/(total(up[7:9,*,0],1))[2],3)  ,linestyle=3,color=100  
oplot,smooth(total(up[10:12,*,0],1)/(total(up[10:12,*,0],1))[3],3)  ,linestyle=4,color=150  
oplot,smooth(total(up[13:15,*,0],1)/(total(up[13:15,*,0],1))[3],3)  ,linestyle=5,color=250  

window,1
loadct,0
plot,total(down[1:3,*,0],1)/(total(down[1:3,*,0],1))[3],color=0,/nodata
loadct,13 
plot,total(down[1:3,*,0],1)      ,color=0         
oplot,total(down[4:6,*,0],1) ,linestyle=2,color=50  
oplot,total(down[7:9,*,0],1) ,linestyle=3,color=100  
oplot,total(down[10:12,*,0],1) ,linestyle=4,color=150  
oplot,total(down[13:15,*,0],1) ,linestyle=5,color=250  


global_pow=fltarr(44,npt)
global_pow2=fltarr(44,npt)
FOR i=0, npt-1 DO BEGIN
    series=reform(vmap_pro[*,i,7,0])
    wave = WAVELET(series,60.,period=period)
    global_pow[*,i]=total(abs(wave)^2,1)/nt

    series=reform(vmap_ret[*,i,7,0])
    wave = WAVELET(series,30.,period=period)
    global_pow2[*,i]=total(abs(wave)^2,1)/nt
END
stop

END