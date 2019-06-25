;PURPOSE - Calculate the non-thermal widths form the Doppler widths


PRO comp_non_thermal,date,nfiles=nfiles

inpath  = 'data/CoMP/wave_tracking_output/'+date+'/'
restore,inpath+'cube_ivw_'+date+'.sav'


IF n_elements(nfiles) EQ 0 THEN nfiles=1
sz=size(cube_w)
nt=sz(3)
IF nfiles GT nt THEN nfiles=nt-1

;Thermal width
;sig_th=sqrt(2*k_b*T_e/m_ion)
;Assuming ion and electron temperatures are the same
;in the CoMP corona, T_e~1.6 MK - peak formation temp
;of Fe XIII

sig_th=21.0 ;21 km/s


;Instrumental width
sig_inst=21. ;21 km/s

cube_nt=sqrt(cube_w[*,*,0:nfiles-1]^2-sig_inst^2-sig_th^2)

save,cube_nt,filename=inpath+'cube_nth_'+date+'.sav'

END