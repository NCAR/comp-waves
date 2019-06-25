PRO mcsimpl_many


outpath='analysis/CoMP/simulations/'
files=find_files('sim_res*.idlsave',outpath)
nfiles=n_elements(files)
restore,files[0]
res_list=list(res)
nt_comb=fltarr(nfiles)
nt_comb[0]=nt
FOR i=1, nfiles-1 DO BEGIN
	restore,files[i]
	res_list.add,res
	nt_comb[i]=nt
ENDFOR

delt_aic=list(res_list[0].aic_line-res_list[0].aic_pg)
FOR i=1,nfiles-1 DO delt_aic.add,res_list[i].aic_line-res_list[i].aic_pg

in=list(where(delt_aic[0] gt 0))
FOR i=1,nfiles-1 DO in.add,where(delt_aic[i] gt 0)

const=list(res_list[0].fit_pg[0])
FOR i=1,nfiles-1 DO const.add,res_list[i].fit_pg[0]

mn_const=fltarr(nfiles,2)
FOR i=0,nfiles-1 DO mn_const[i,*]=(moment(alog((const[i])[in[i]])))[0:1]

slope=list(res_list[0].fit_pg[1])
FOR i=1,nfiles-1 DO slope.add,res_list[i].fit_pg[1]

mn_slope=fltarr(nfiles,2)
FOR i=0,nfiles-1 DO mn_slope[i,*]=(moment((slope[i])[in[i]]))[0:1]

noise=list(res_list[0].fit_pg[2])
FOR i=1,nfiles-1 DO noise.add,res_list[i].fit_pg[2]

mn_noise=fltarr(nfiles,2)
FOR i=0,nfiles-1 DO mn_noise[i,*]=(moment((noise[i])[in[i]]))[0:1]

amp=list(res_list[0].fit_pg[3])
FOR i=1,nfiles-1 DO amp.add,res_list[i].fit_pg[3]

mn_amp=fltarr(nfiles,2)
FOR i=0,nfiles-1 DO mn_amp[i,*]=(moment((amp[i])[in[i]]))[0:1]


cent_freq=list(res_list[0].fit_pg[4])
FOR i=1,nfiles-1 DO cent_freq.add,res_list[i].fit_pg[4]

mn_cent_freq=fltarr(nfiles,2)
FOR i=0,nfiles-1 DO mn_cent_freq[i,*]=(moment((cent_freq[i])[in[i]]))[0:1]

freq_wid=list(res_list[0].fit_pg[5])
FOR i=1,nfiles-1 DO freq_wid.add,res_list[i].fit_pg[5]

mn_freq_wid=fltarr(nfiles,2)
FOR i=0,nfiles-1 DO mn_freq_wid[i,*]=(moment(abs((freq_wid[i])[in[i]])))[0:1]

p_main=errorplot(nt_comb,mn_const[*,0],sqrt(mn_const[*,1]),'2b+',layout=[1,6,1],ytickint=0.5,$
	             xtickform="(A1)", pos=[0.15,0.84,0.85,0.98],ytitle='log(A)',sym_thick=2,$
	             xthick=2,ythick=2,ymin=0,/xlog,xr=[100,3000],xst=1,/buffer)
line=plot([100,3000],alog(p[0])*[1,1],'--1r',overplot=p_main)
p_2=errorplot(nt_comb,mn_slope[*,0],sqrt(mn_slope[*,1]),'2b+',layout=[1,6,2],/current,$
			ytickint=0.1,xtickform="(A1)",pos=[0.15,0.69,0.85,0.83],yran=[-1.35,-1.05],sym_thick=2,$
	             xthick=2,ythick=2,ymin=0,/xlog,xr=[100,3000],xst=1)
line=plot([100,3000],(p[1])*[1,1],'--1r',overplot=p_2)
text=text(0.065,0.74,'$\alpha$',orient=90)
p_3=errorplot(nt_comb,mn_noise[*,0]*1e3,sqrt(mn_noise[*,1])*1e3,'2b+',layout=[1,6,3],/current, $
	      xtickform="(A1)",pos=[0.15,0.54,0.85,0.68],yran=[0.4,1.7],ytickval=[1.0,1.5],sym_thick=2,$
	             xthick=2,ythick=2,ymin=0,/xlog,xr=[100,3000],xst=1)
line=plot([100,3000],(p[2])*[1,1]*1e3,'--1r',overplot=p_3)
text=text(0.065,0.59,'B',orient=90)
p_4=errorplot(nt_comb,mn_amp[*,0]*1e3,sqrt(mn_amp[*,1])*1e3,'2b+',layout=[1,6,4],/current, $
	        ytickint=1.,xtickform="(A1)",pos=[0.15,0.39,0.85,0.53],sym_thick=2,$
	             xthick=2,ythick=2,ymin=0,/xlog,xr=[100,3000],xst=1)
line=plot([100,3000],(p[3])*[1,1]*1e3,'--1r',overplot=p_4)
text=text(0.065,0.44,'C',orient=90)
p_6=errorplot(nt_comb,mn_cent_freq[*,0],sqrt(mn_cent_freq[*,1]),'2b+',layout=[1,6,5],/current,ytickint=0.1,$
				xtickform="(A1)",pos=[0.15,0.24,0.85,0.38],sym_thick=2,$
	             xthick=2,ythick=2,ymin=0,/xlog,xr=[100,3000],xst=1)
tex=text(0.065,0.29,'D',orient=90)
line=plot([100,3000],(p[4])*[1,1],'--1r',overplot=p_6)
p_5=errorplot(nt_comb,mn_freq_wid[*,0],sqrt(mn_freq_wid[*,1]),'2b+',layout=[1,6,6],/current,ytickint=0.1, $
				pos=[0.15,0.09,0.85,0.23],xtitle='N',sym_thick=2,$
	             xthick=2,ythick=2,ymin=0,/xlog,xr=[100,3000],xst=1)
line=plot([100,3000],(p[5])*[1,1],'--1r',overplot=p_5)
text=text(0.065,0.14,'E',orient=90)
stop
p_main.Save,outpath+'sim_res_plot_all.png'
stop


END