;PURPOSE: 2-D windowed padded FFT
;         Applies two corrections to power for padding and windowing
;         Returns all frequencies and positive wavenumbers
;
;INPUTS: inarr - 2d array to be FFT'd
;        dt - cadence
;        dx - resolution
;        
;OPTIONAL INPUTS - ybas - controls elements in zero padding in space, default is 0
;                - tbas - controls elements in zero padding in time, default is 0
;
;OUTPUTS freq_n - frequency array
;        k_n - wavenumber array
;
;
;Created R J Morton MAY 2015
;
;
;

FUNCTION twod_fft, inarr,dt,dx,freq_n=freq_n,k_n=k_n,ybas=ybas,tbas=tbas,han=han

sz=size(inarr)
ny=sz(1) & nt=sz(2)

IF n_elements(han) eq 0 THEN hn=hanning(ny,nt,alpha=0.54) ELSE BEGIN
    ; define temporal apodization
    apodt = fltarr(nt)+1
    apod=0.1
    apodrimt = nt*apod
    apodt[0] = (sin(!pi/2.*findgen(apodrimt)/apodrimt))^2
    apodt = apodt*shift(rotate(apodt,2),1) 

    ; define spatial apodization
    apody = fltarr(ny)+1
    apodrimy = ny*apod
    apody[0] = (sin(!pi/2.*findgen(apodrimy)/apodrimy))^2
    apody = apody*shift(rotate(apody,2),1) 

    hn=transpose(rebin(apodt,nt,ny))*rebin(apody,ny,nt)

ENDELSE

IF n_elements(ybas) EQ 0 THEN ybas=0
IF n_elements(tbas) EQ 0 THEN tbas=0

IF ybas GT 0 THEN base2 = FIX( ALOG(ny) / ALOG(2) ) + ybas   ; power of 2 nearest to N
IF tbas GT 0 THEN base2t = FIX( ALOG(nt) / ALOG(2) ) + tbas   ; power of 2 nearest to N

IF ybas GT 0 THEN yel=2L^(base2) ELSE yel=ny
IF tbas GT 0 THEN tel=2L^(base2t) ELSE tel=nt

s_tot=n_elements(inarr)
cpg=total(hn)/s_tot ;cpg - coherent power gain - correction factor for apodisation

		
ts=fltarr(yel,tel)
ts[0:ny-1,0:nt-1]=inarr*hn

n_tot=n_elements(ts)
n_len=size(ts)
fac=(1.*n_tot/s_tot)^2 ; Correction factor required for padding
pow=fft(ts)*sqrt(fac/cpg^2)


half=(n_len[2]-1)/2
X = FINDGEN(half) + 1
IF (n_len[2] MOD 2) EQ 0 THEN $
  freq_n = [0.0, X, n_len[2]/2, -n_len[2]/2 + X]/(n_len[2]*dt) $
ELSE $
  freq_n = [0.0, X, -(n_len[2]/2 + 1) + X]/(n_len[2]*dt)

half=(n_len[1]-1)/2
X = FINDGEN(half) + 1
IF (n_len[1] MOD 2) EQ 0 THEN $
  k_n = [0.0, X]/(n_len[1]*dx) $
ELSE $
  k_n = [0.0, X]/(n_len[1]*dx)
;k_n=findgen(n_len[1]/2)/(n_len[1]*dx)


return,pow

END