pro wave_mc

nreps=1000
dt=30
nt=164

;For setup
ser=randomn(systime,nt)
pow=abs(wavelet(ser,dt,period=period,/pad,scale=scale,j=j1,$
                      s0=s0,coi=coi,dj=dj) )^2


nscale=n_elements(scale)


stop
;Define mask for removing COI data
mask=fltarr(nt,nscale,/nozero)
mask[*]=1.
FOR i=0,nt-1 DO BEGIN
   ht=where(period ge coi[i])
   mask(i,ht)=0.
ENDFOR

;95% point estimate confidence level
sig= -alog(0.05)

scl=0.
cnt_pix=0.
wvpacks=0.

;window,0

;Increase series length by 3 orders (2^10)
base2 = FIX(ALOG(nt)/ALOG(2) + 0.4999) + 10
nt2=2L^base2

;Set power spectra
freq = (findgen(nt2/2)+1)/nt2/dt
p=[8.5e-6,-1.2,0.001,0.0034,-5.57,0.32]
pow=p[0]*freq^p[1]

;Generate long time-series
x = TS_GEN(nt2, dt=dt, freq=freq, pow=pow,time=time,/quiet)

; define temporal apodization
apodt = fltarr(nt)+1
apod=0.1
apodrimt = nt*apod
apodt[0] = (sin(!pi/2.*findgen(apodrimt)/apodrimt))^2
apodt = apodt*shift(rotate(apodt,2),1) 

;create random array indexs
index=fix(randomu(1,nreps)*(nt2-nt-1),type=3)
coeff_fits=fltarr(nreps,2)
perc=list(0)


FOR i=0,nreps-1 DO BEGIN

	ser=x[index[i]:index[i]+nt-1]
    input=(ser-mean(ser))

    pow_coeff=mle_fit_psd(input*apodt,dt,freqo,pergm)
    coeff_fits[i,*]=pow_coeff

    back_fit_per=double(10.^(pow_coeff[0]+pow_coeff[1]*ALOG10(1./period)))
    back_fit_per=back_fit_per/dt/2. ;Correct Vaughan normalisation for comparison to wavelet
    
    pow=abs(wavelet(input,dt,period=period,/pad,scale=scale,j=j1,$
                      s0=s0,coi=coi,dj=dj) )^2

    pow=pow*mask
    norma=pow/transpose(rebin(back_fit_per,nscale,nt))
    ;IF n_elements(pow_dt) EQ 0 THEN pow_dt=norma $
    ;ELSE pow_dt=[[[temporary(pow_dt)]],[[norma]]]

    lb=label_region(norma gt sig)
     
    perc.add,n_elements(where(lb gt 0))*1./nt/nscale 

    ;go through blobs
    FOR j=1, max(lb) DO BEGIN
    	in=where(lb eq j) ; find blob 
    	in2=in/nt         ; find which scale each blob row is on
    	in3=in mod nt
        pdf=histogram(in2,locations=xbin) ; find number of pixels at each scale
        scl=[temporary(scl),xbin]
        cnt_pix=[temporary(cnt_pix),pdf]
        wvpacks=[temporary(wvpacks),dt/period[xbin]*pdf]
        ;window,1
    	;plot,
    	
    ENDFOR
    ;stop
    counter,i,nreps,'Calculating',/percent
ENDFOR

;mn=fltarr(max(scl)+1)
;ninefive=fltarr(nscale)
;FOR i=0, nscale-1 DO ninefive[i]=percentiles(pow_dt[where(mask[*,i] gt 0),i,*],value=0.95)
;    in=where(scl eq i)

;    edf,cnt_pix[in],x,y

;    lvl=where(y lt 0.95)
;    ninefive[i]=x[lvl[-1]+1]
;    mn[i]=mean(cnt_pix[in])
;ENDFOR

stop
END