PRO CH_coher

date='20120327'

inpath  = 'analysis/CoMP/wave_tracking_output/'+date+'/'
restore,inpath+'boxes_for_average.sav'

cad=30. ;Temporal cadence

sz=size(box)

box_sep=fltarr(sz(1),sz(2),sz(3),2)

NoYs = (sz(2) - 1)/2 + 1 ;Will run from 0 freq to middle - excluding Nyquist if even
NoTs = (sz(3) - 1)/2 + 1

IF sz(2) mod 2 EQ 0 THEN exY=1 ELSE exY=0
IF sz(3) mod 2 EQ 0 THEN exT=1 ELSE exT=0

;separate inward and outward
FOR i=0,sz(1)-1 DO BEGIN
	ft=fft(reform(box[i,*,*]))
    outward=ft
    outward[0:NoYs,0:NoTs]=0. & outward[NoYs+1+exY:-1,NoTs+1+exT:-1]=0.
    inward=ft
    inward[0:NoYs,NoTs+1+exT:-1]=0. & inward[NoYs+1+exY:-1,0:NoTs]=0.
    box_sep[i,*,*,0]=fft(outward,/inverse)
    box_sep[i,*,*,1]=fft(inward,/inverse)
ENDFOR


length=60
overlap=0.3
range=7

;test to get right dimensions
ps_welch,reform(box[0,0,*]),length=length,ft_ser=ft_ser,overlap=overlap,/verb

ftsz=size(ft_ser)
ft_out_len=ftsz(1)/2
ft_outs=complexarr(sz(1),sz(2),ft_out_len,ftsz(2)) ;size is nx, ny, nt/2, no series


FOR i=0,sz(1)-1 DO FOR j=0,sz(2)-1 DO BEGIN

ps_welch,reform(box_sep[i,j,*,1]),length=length,ft_ser=ft_ser,overlap=overlap
IF overlap gt 0 THEN ft_outs[i,j,*,*]=ft_ser[0:ft_out_len-1,*] ELSE ft_outs[i,j,*]=ft_ser[0:ft_out_len-1] 

ENDFOR


;coher=fltarr(ft_out_len,range,sz(2)*(sz(1)-range+1))


temp={spec_ot:complexarr(ft_out_len,ftsz(2),range), spec_oo:fltarr(ft_out_len,ftsz(2),range), $
      spec_tt:fltarr(ft_out_len,ftsz(2),range)}


FOR i=0,sz(2)-1 DO BEGIN ; loop over x
  FOR j=range/2,sz(1)-range/2-1 DO BEGIN ;loop over y avoiding edges
      ref=reform(ft_outs[j,i,*,*])
  	FOR k=-range/2,range/2 DO BEGIN ;compare series at (i,j) to  nearest neighbours
		
		comp=reform(ft_outs[j+k,i,*,*])
		temp.spec_ot[*,*,k+range/2]=ref*conj(comp)
		temp.spec_oo[*,*,k+range/2]=real_part(ref*conj(ref))
		temp.spec_tt[*,*,k+range/2]=real_part(comp*conj(comp))

        ;coher[*,k+range/2,h]=abs(mean(ref*conj(comp),dim=2))^2/real_part(mean(ref*conj(ref),dim=2))/$
	    ;          real_part(mean(comp*conj(comp),dim=2))
    ENDFOR
    
    IF n_elements(res) EQ 0 THEN res=temp ELSE res=[temporary(res),temp]

  ENDFOR
ENDFOR

IF overlap GT 0 THEN BEGIN
	xspec_12=mean(mean(res.spec_ot,dim=4),dim=2)
	xspec_11=mean(mean(res.spec_oo,dim=4),dim=2)
	xspec_22=mean(mean(res.spec_tt,dim=4),dim=2)
	coher=abs(xspec_12)/sqrt(xspec_11*xspec_22)    ;coherence

	;Bootstrap estimates for SD
	nelm_dat=(sz(1)-range+1)*sz(2)
	num_bs=500 ;nelm_dat
	mean_coher_bs=fltarr(ft_out_len,range,num_bs)
	bs_coher_details=fltarr(ft_out_len,range,2)
	;test_coher_bs=complexarr(ft_out_len,range,num_bs,3)
	index_bs=round(randomu(1,nelm_dat,num_bs)*(nelm_dat-1))
    
    FOR k=0,num_bs-1 DO BEGIN 
        
    	xspec_12=mean(mean((reform(res.spec_ot))[*,*,*,index_bs[*,k]],dim=4),dim=2)
		xspec_11=mean(mean((reform(res.spec_oo))[*,*,*,index_bs[*,k]],dim=4),dim=2)
		xspec_22=mean(mean((reform(res.spec_tt))[*,*,*,index_bs[*,k]],dim=4),dim=2)

		mean_coher_bs[*,*,k]=abs(xspec_12)/sqrt(xspec_11*xspec_22)    ;coherence
		
	ENDFOR

	resolve_routine, 'gauss_cdf',/is_function
    FOR i=0,ft_out_len-1 DO FOR j=0,range-1 DO BEGIN
    	IF j ne range/2 THEN BEGIN
    	    mom=moment(mean_coher_bs[i,j,*],dim=3)
    	    ksone,(mean_coher_bs[i,j,*]-mom[0])/sqrt(mom[1]),'gauss_cdf',d,prob
    		bs_coher_details[i,j,*]=[mom[1]^0.5,d]
    	ENDIF

    ENDFOR

ENDIF ELSE BEGIN

	spec_ot=reform(res.spec_ot[*,0,*,*])
    spec_oo=reform(res.spec_oo[*,0,*,*])
    spec_tt=reform(res.spec_tt[*,0,*,*])
	
	xspec_12=(mean(spec_ot,dim=3))
	xspec_11=(mean(spec_oo,dim=3))
	xspec_22=(mean(spec_tt,dim=3))
	coher=abs(xspec_12)/sqrt(xspec_11*xspec_22)    ;coherence

	;Bootstrap estimates for SD
	nelm_dat=(sz(1)-range+1)*sz(2)
	num_bs=500 ;nelm_dat
	mean_coher_bs=fltarr(ft_out_len,range,num_bs)
	bs_coher_details=fltarr(ft_out_len,range,2)
	;test_coher_bs=complexarr(ft_out_len,range,num_bs,3)
	index_bs=round(randomu(1,nelm_dat,num_bs)*(nelm_dat-1))
    
    spec_ot=reform(res.spec_ot[*,0,*,*])
    spec_oo=reform(res.spec_oo[*,0,*,*])
    spec_tt=reform(res.spec_tt[*,0,*,*])

    FOR k=0,num_bs-1 DO BEGIN 
    ;    CASE k OF 
    ;      0: index_jk=indgen(num_bs-1)+1
    ;      num_bs-1: index_jk=indgen(num_bs-1)
    ;      ELSE: index_jk=[indgen(k),indgen(num_bs-1-k)+k+1]
	;	ENDCASE
        
    	xspec_12=(mean(spec_ot[*,*,index_bs[*,k]],dim=3))
		xspec_11=(mean(spec_oo[*,*,index_bs[*,k]],dim=3))
		xspec_22=(mean(spec_tt[*,*,index_bs[*,k]],dim=3))

	;	xspec_12=(mean(spec_ot[*,*,index_jk],dim=3))
	;	xspec_11=(mean(spec_oo[*,*,index_jk],dim=3))
;		xspec_22=(mean(spec_tt[*,*,index_jk],dim=3))
;		test_coher_bs[*,*,k,0:2]=[xspec_12,xspec_11,xspec_22]
		mean_coher_bs[*,*,k]=abs(xspec_12)/sqrt(xspec_11*xspec_22)    ;coherence
		
	ENDFOR

	resolve_routine, 'gauss_cdf',/is_function
    FOR i=0,ft_out_len-1 DO FOR j=0,range-1 DO BEGIN
    	IF j ne range/2 THEN BEGIN
    	    mom=moment(mean_coher_bs[i,j,*],dim=3)
    	    ksone,(mean_coher_bs[i,j,*]-mom[0])/sqrt(mom[1]),'gauss_cdf',d,prob
    		bs_coher_details[i,j,*]=[mom[1]^0.5,d]
    	ENDIF

    ENDFOR

ENDELSE


freq=findgen(ft_out_len)/ft_out_len/2./cad

;Plots of coherence for different separations
pl=plot(freq*1e3,coher[*,range/2],'Tu',yr=[0,1.1],xr=[0,max(freq)]*1e3,xtitle='Frequency (mHz)',ytitle='Coherence',name='Zero')
pl2=plot(freq*1e3,coher[*,range/2-1],'r',overplot=pl,name='3.2 Mm l')
pl3=plot(freq*1e3,coher[*,range/2-2],'b',overplot=pl,name='6.5 Mm l')
pl4=plot(freq*1e3,coher[*,range/2-3],'g',overplot=pl,name='9.7 Mm l')
pl5=plot(freq*1e3,coher[*,range/2+1],'r--',overplot=pl,name='3.2 Mm r')
pl6=plot(freq*1e3,coher[*,range/2+2],'b--',overplot=pl,name='6.5 Mm r')
pl7=plot(freq*1e3,coher[*,range/2+3],'g--',overplot=pl,name='9.7 Mm r')
leg=legend(target=[pl,pl2,pl3,pl4,pl5,pl6,pl7],pos=[0.9,0.87])

;Find exponential fits for correlation length

x=[-reverse(findgen(range/2)-1.),findgen(range/2+1)]*3.2263258 ; Mm
vals={fit_results, params:fltarr(2), param_errs:fltarr(2), dof:1.,chisq:1.}
vals=replicate(vals,ft_out_len)

;weights=replicate(1.,range)
weights=1./bs_coher_details[*,*,0]^2
A=[0.5,0.2]
FOR i=0,ft_out_len-1 DO BEGIN
  ys=reform(coher[i,*])
  errs=reform(bs_coher_details[i,*,0])
  errs[range/2]=errs[range/2-1]
 ; help,ys,errors
  ress=mpfitfun('mycauchy',x,ys,errs,A,/quiet,perror=perror,dof=dof,bestnorm=bestnorm)
  ;A=[0.5,0.2]
  ;ress=curvefit(x,reform(coher[i,*]),reform(weights[i,*]),A,sigma,function_name='mycauchy_pro')
  
     vals[i].params=ress
     vals[i].param_errs=perror
     vals[i].dof=dof
     vals[i].chisq=bestnorm
 ; ENDIF ELSE vals=[vals,{params:ress, param_errs:perror, dof:dof,chisq:bestnorm}]
ENDFOR

cal_err=abs(vals[*].param_errs[1])/abs(vals[*].params[1])^2
cor_len_pl=errorplot(freq*1e3,1./abs(vals[*].params[1]),cal_err,xtitle='Frequency (mHz)',ytitle='Correlation length (Mm)',xr=[0,max(freq)]*1e3)
res_pl=plot([freq(0),freq(-1)]*1e3,[1,1]*3.2263258,overplot=cor_len_pl,linestyle=2,color='r',thick=2)


stop

END