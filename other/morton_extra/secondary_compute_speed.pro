FUNCTION secondary_compute_speed,main,comp,freq

s_tot=n_elements(main)
h=hanning(s_tot,alpha=0.54)
cpg=total(h)/s_tot

base2t = FIX( ALOG(s_tot) / ALOG(2) ) + 4  ; power of 2 nearest to N
pad_main=fltarr(2L^(base2t))
pad_comp=fltarr(2L^(base2t))
pad_main[0:s_tot-1]=main*h
pad_comp[0:s_tot-1]=comp*h
n_tot=n_elements(pad_main)

fac=(1.*n_tot/s_tot)^2

freq=findgen(n_tot)/n_tot/30.
freq0=0.003
fwidth = 0.001
filter = exp( -(freq-freq0)^2/fwidth^2 )
filter(0) = 0.                      ; set dc to zero
filter(where (freq lt .001)) = 0.   ; set low frequencies to zero
filter=filter/total(filter)


cp=fft(pad_main)*conj(fft(pad_comp))
wcp=filter*abs(fft(pad_main)*conj(fft(pad_comp)))*fac/cpg^2
phi=atan(imaginary(cp)/real_part(cp))
a=max(wcp[1:n_tot/2],loc)
phase=phi(loc)
lag=phase/2./!PI/freq(loc)

return,[lag,2.*sqrt(a)]

END