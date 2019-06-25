
;creates random distribution of kink waves from known SDO measurements
PRO rand_kw

restore,'analysis/sdo/2012/03/27/20120327_1600-2000_east_qs_50Mm_mean_nuwt_vel_amp.sav'
nuwt_va_qs=nuwt_vel_amp

in=where(nuwt_va_qs.peak_freq lt 0.03) ; remove very high frequency values

freqs=nuwt_va_qs.peak_freq[in]
vamps=nuwt_va_qs.peak_vel_amp[in]


;work in log space where distributions are approximately gaussian
freq_mom=moment(alog10(freqs))
vamps_mom=moment(alog10(vamps))

cov=correlate(alog10(freqs),alog10(vamps),/covariance)

covar=[[freq_mom[1],cov],[cov,vamps_mom[1]]]
stop
out=mrandomn(1,covar,700)
out2=out+transpose(rebin([freq_mom[0],vamps_mom[0]],2,700))

out3=10^out2

;create a time-series
ser=fltarr(5000)

;random start time
t0=round(randomu(1,239)*4999)

;length of signal in periods
length=nuwt_VA_QS.peaK_duration[in]*nuwt_va_qs.peak_freq[in]
length_osc=length[round(randomu(2,239)*(n_elements(freqs)-1))]

;this doesn't work
FOR i=0,238 DO BEGIN

	wave=out3[i,1]*sin(findgen(round(length_osc[i]/out3[i,0]))*2.*!pi*out3[i,0])
    ln=n_elements(wave)
    endp=t0[i]+ln
    over_ln=5000-endp
    IF over_ln lt 0 then ser[t0[i]:-1]+=wave[0:over_ln] ELSE ser[t0[i]:t0[i]+ln-1]+=wave

ENDFOR


stop

END


