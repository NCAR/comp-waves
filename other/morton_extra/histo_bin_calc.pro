PRO histo_bin_calc,x,minbin,maxbin

N=n_elements(x)


eval=fltarr(maxbin-minbin+1,4)
FOR alp=0,3 DO BEGIN
  
    alpha=10.^(alp-2)
    
	FOR ii=minbin,maxbin-1 DO BEGIN
	    
		lhood=0
		xnew=x
		bs=ii
        plothist,xnew,bin=bs,xhist,yhist,/noplot
        ;FOR jj=0,N-1 DO BEGIN
		;     xnew=shift(xnew,1)
	    ;     bin_ind=floor( (xnew[0]-min(xnew))/ii)
	         ;print,(yhist[bin_ind]+alpha-1)/ii/(total(yhist+alpha)-1) 
	         ;pause
	    ;     lhood=[lhood,yhist(bin_ind)*alog( (yhist[bin_ind]+alpha-1.)/ii/(total(yhist+alpha)-1.) )]
        ;ENDFOR
        ; print,lhood
        ; pause
         ;help,lhood,ii
         eval[ii-minbin,alp]= total(yhist*alog( (yhist+alpha-1.)/ii/(total(yhist+alpha)-1.) ));total(lhood[1:N])
        
	ENDFOR


ENDFOR
stop
END