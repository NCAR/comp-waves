PRO quick_speed_compute,date,debug=debug,name_addon=name_addon

 
file='data/comp/wave_tracking_output/'+date+'/'
outpath='data/comp/wave_tracking_output/'+date+'/'

restore,file+'cube_ivw_'+date+'.sav',/verbose
restore,file+'wave_angle_'+date+'_3_3.5mHz.sav',/verbose

sz=size(cube_v)
nx=sz(1) & ny=sz(2) & nt=sz(3)

lower_r = 231.
upper_r = 270.
nwlst = strcompress(string(3),/rem)


maxoffset=0.
if n_elements(name_addon) eq 0 then begin
 name_addon = ''
endif

if keyword_set(wr_speeds) then wr_speeds = 1 else wr_speeds = 0

;diff=cube_v-smooth(cube_v,[1,1,3],/edge_truncate)
;diff=sqrt(mean(diff^2,dim=3))
;in=where(diff gt 0.5)

;FOR i=0,nt-1 DO BEGIN 
;  dum=cube_v[*,*,i]
;  dum[in]=0.
;  cube_v[*,*,i]=dum
;ENDFOR

;cube_v=cube_v[*,*,0:50]
wave_angle=median(wave_angle,3)



space_time_run, date, cube_v, wave_angle, index, $
                 outpath, nwlst, '3.5', name_addon, debug=debug,wr_speeds=wr_speeds,$
                 npt_init=25,r_upper=upper_r,r_inner=lower_r


END