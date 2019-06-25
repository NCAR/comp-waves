PRO noise_spectra,date,file_ext=file_ext,err_fit=err_fit,doap=doap


inpath  = 'analysis/CoMP/wave_tracking_output/'+date+'/'
restore,inpath+file_ext+'.sav'

dt=30.

me=mean(err_dv,dimension=3)
sz=size(err_dv)
nx=sz(1)
ny=sz(2)
nt=sz(3)


ferr=findgen(nt/2)/(nt*dt)

gnoise=randomn(systime_seed,nx,ny,nt)

noise=fltarr(nx,ny,nt)

FOR i=0,nt-1 DO noise[*,*,i]=gnoise[*,*,i]*me

err_spec=complexarr(nx,ny,nt)

IF keyword_set(doap) THEN BEGIN
   hn=hanning(nt)
   cpg=total(hn)/nt
ENDIF ELSE BEGIN
   hn=1.
   cpg=1.
ENDELSE

for i=0,sz(1)-1 do begin                 
for j=0,sz(2)-1 do begin
  in=noise[i,j,*]-mean(noise[i,j,*])                 
  err_spec[i,j,*]=fft(in*hn,-1)      
endfor
endfor

err_tot_spec=total(total(abs(err_spec[*,*,1:nt/2-1])^2/cpg^2,2),1)/(nx*ny)
res=linfit(findgen(nt/2-1),err_tot_spec)
err_fit=res[0]+res[1]*findgen(nt/2-2)



END