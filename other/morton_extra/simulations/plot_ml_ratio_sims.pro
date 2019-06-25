PRO plot_ml_ratio_sims

inpath='analysis/CoMP/simulations/'

restore,inpath+'sim_res_ratio_N10_L160.idlsave'
res_n10_l160=res
mn_res_n10_l160=moment(res_n10_l160.fit,dim=2)
restore,inpath+'sim_res_ratio_N10_L320.idlsave'
res_n10_l320=res
mn_res_n10_l320=moment(res_n10_l320.fit,dim=2)
restore,inpath+'sim_res_ratio_N10_L640.idlsave'
res_n10_l640=res
mn_res_n10_l640=moment(res_n10_l640.fit,dim=2)
restore,inpath+'sim_res_ratio_N10_L960.idlsave'
res_n10_l960=res
mn_res_n10_l960=moment(res_n10_l960.fit,dim=2)
restore,inpath+'sim_res_ratio_N10_L1280.idlsave'
res_n10_l1280=res
mn_res_n10_l1280=moment(res_n10_l1280.fit,dim=2)
restore,inpath+'sim_res_ratio_N10_L1600.idlsave'
res_n10_l1600=res
mn_res_n10_l1600=moment(res_n10_l1600.fit,dim=2)

meanvals=[mn_res_n10_l160[0,0],mn_res_n10_l320[0,0],mn_res_n10_l640[0,0],mn_res_n10_l960[0,0],mn_res_n10_l1280[0,0],mn_res_n10_l1600[0,0]]
meanvals_var=[mn_res_n10_l160[0,1],mn_res_n10_l320[0,1],mn_res_n10_l640[0,1],mn_res_n10_l960[0,1],mn_res_n10_l1280[0,1],mn_res_n10_l1600[0,1]]
xvals=[1.,2.,4.,6.,8.,10]*160
pl=errorplot(xvals,(meanvals-1.0)*100.,100.*meanvals_var^0.5/sqrt(5000.),xrange=[0,xvals[-1]+50],$
				ytitle='$b_{p_0} (%)$',/buffer,line=6,symbol='circle',errorbar_capsize=0,$
				layout=[1,2,1],pos=[0.15,0.49,0.85,0.83], xtickform="(A1)")
pl2=plot([0,xvals[-1]+50],[0.,0.],'--',/overplot)

meanvals2=[mn_res_n10_l160[1,0],mn_res_n10_l320[1,0],mn_res_n10_l640[1,0],mn_res_n10_l960[1,0],mn_res_n10_l1280[1,0],mn_res_n10_l1600[1,0]]
meanvals_var2=[mn_res_n10_l160[1,1],mn_res_n10_l320[1,1],mn_res_n10_l640[1,1],mn_res_n10_l960[1,1],mn_res_n10_l1280[1,1],mn_res_n10_l1600[1,1]]
pl3=errorplot(xvals,(meanvals2-100.0),meanvals_var2^0.5/sqrt(5000.),xrange=[0,xvals[-1]+50],$
		     xtitle='Length of series',ytitle='$b_{p_1} (%)$',/current,line=6,symbol='circle',errorbar_capsize=0,$
				layout=[1,2,2],pos=[0.15,0.11,0.85,0.45])
pl4=plot([0,xvals[-1]+50],[0.,0.],'--',/overplot)
pl.Save,inpath+'plot1.png'

restore,inpath+'sim_res_ratio_N1_L160.idlsave'
res_n1_l160=res[where(res.fit[0] ne 0 )]
mn_res_n1_l160=moment(res_n1_l160.fit,dim=2)
restore,inpath+'sim_res_ratio_N5_L160.idlsave'
res_n5_l160=res[where(res.fit[0] ne 0 )]
mn_res_n5_l160=moment(res_n5_l160.fit,dim=2)
restore,inpath+'sim_res_ratio_N20_L160.idlsave'
res_n20_l160=res[where(res.fit[0] ne 0 )]
mn_res_n20_l160=moment(res_n20_l160.fit,dim=2)
restore,inpath+'sim_res_ratio_N30_L160.idlsave'
res_n30_l160=res[where(res.fit[0] ne 0 )]
mn_res_n30_l160=moment(res_n30_l160.fit,dim=2)
restore,inpath+'sim_res_ratio_N50_L160.idlsave'
res_n50_l160=res[where(res.fit[0] ne 0 )]
mn_res_n50_l160=moment(res_n50_l160.fit,dim=2)

meanvals3=[mn_res_n1_l160[0,0],mn_res_n5_l160[0,0],mn_res_n10_l160[0,0],mn_res_n20_l160[0,0],mn_res_n30_l160[0,0],mn_res_n50_l160[0,0]]
meanvals_var3=[mn_res_n1_l160[0,1],mn_res_n5_l160[0,1],mn_res_n10_l160[0,1],mn_res_n20_l160[0,1],mn_res_n30_l160[0,1],mn_res_n50_l160[0,1]]
xvals=[1.,2.,10.,20.,30.,50.]
pl=errorplot(xvals,(meanvals3-1.0)*100.,100.*meanvals_var3^0.5/sqrt(5000.),xrange=[0,xvals[-1]+10],yrange=[-1,10],$
				ytitle='$b_{p_0} (%)$',/buffer,line=6,symbol='circle',errorbar_capsize=0,$
				layout=[1,2,1],pos=[0.15,0.49,0.85,0.83], xtickform="(A1)")
pl2=plot([0,xvals[-1]+10],[0.,0.],'--',/overplot)


meanvals4=[mn_res_n1_l160[1,0],mn_res_n5_l160[1,0],mn_res_n10_l160[1,0],mn_res_n20_l160[1,0],mn_res_n30_l160[1,0],mn_res_n50_l160[1,0]]
meanvals_var4=[mn_res_n1_l160[1,1],mn_res_n5_l160[1,1],mn_res_n10_l160[1,1],mn_res_n20_l160[1,1],mn_res_n30_l160[1,1],mn_res_n50_l160[1,1]]
pl3=errorplot(xvals,(meanvals4-100.0),meanvals_var4^0.5/sqrt(5000.),xrange=[0,xvals[-1]+10],$
				xtitle='Number of series averaged together',ytitle='$b_{p_1} (%)$',/buffer,line=6,symbol='circle',$
				errorbar_capsize=0,layout=[1,2,2],pos=[0.15,0.11,0.85,0.45],/current)
pl4=plot([0,xvals[-1]+10],[0.,0.],'--',/overplot)
pl.Save,inpath+'plot2.png'

stop

END