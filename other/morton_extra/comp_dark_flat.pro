
;Estimate dark/flat noise levels of CoMP
;

PRO comp_dark_flat

files='data/comp/2012/03/27/flats_darks/'

dark_f=files+'20120327.084826.fts'
flat_f=files+'20120327.084943.fts'

files_no=35

im=readfits(dark_f,ext=1,/silent)
sz=size(im)

darks=fltarr(sz(1),sz(2),files_no-1)

;Read in dark files
FOR i=1,files_no-1 DO darks[*,*,i-1]=readfits(dark_f,ext=i,/silent)

diff=darks-shift(darks,[0,0,1])
window,0
plothist,diff,xhist,yhist,xrange=[-50,50],xtitle='dark frame2frame variation'
res=gaussfit(xhist,yhist,coeff,nterms=3)
oplot,xhist,res
print,'sigma',coeff[2]


window,1
;Flats more complicated
flats=fltarr(sz(1),sz(2),files_no-1)
wvleng=strarr(files_no)
sigma=fltarr(files_no,4)

;Read in flat files
FOR i=1,files_no-1 DO BEGIN
   flats[*,*,i-1]=readfits(flat_f,hd,ext=i,/silent)
   wvleng[i]=hd[7]
   
   
ENDFOR
diff=flats[0:610,440:1023,*]-shift(flats[0:610,440:1023,*],[0,0,1])
FOR i=0,files_no-1 DO BEGIN
    in=where(flats[0:610,440:1023,i-1] gt 3000)
    dum=diff[*,*,i-1]
    plothist,dum[in],xhist,yhist,xrange=[-500,500],xtitle='flat frame2frame variation',/noplot
    res=gaussfit(xhist,yhist,coeff,nterms=3)
    ;oplot,xhist,res
    sigma[i,0]=coeff[2]
    sigma[i,2]=median(flats[0:610,440:1023,i-1])
    
ENDFOR

diff=flats[420:1023,0:620,*]-shift(flats[420:1023,0:620,*],[0,0,1])
FOR i=0,files_no-1 DO BEGIN
    in=where(flats[420:1023,0:620,i-1] gt 3000)
    dum=diff[*,*,i-1]
    plothist,dum[in],xhist,yhist,xrange=[-500,500],xtitle='flat frame2frame variation',/noplot
    res=gaussfit(xhist,yhist,coeff,nterms=3)
    ;oplot,xhist,res
    sigma[i,1]=coeff[2]
    sigma[i,3]=median(flats[420:1023,0:620,i-1])
    
ENDFOR
print,median(sigma[*,0]),median(sigma[*,1])
print,sqrt(sigma[*,2:3])

END