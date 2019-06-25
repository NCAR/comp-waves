PRO comp_despike,cube,index,despike


nx=index.NAXIS1
ny=index.NAXIS2
nt=index.NFRAMES

;apply mask
cube=cube*rebin(index.mask,nx,ny,nt)

x=findgen(nx)
y=x

x=rebin(x,nx,ny)
y=transpose(rebin(y,nx,ny))

r=sqrt((x-nx/2)^2+(y-ny/2)^2)

despike=fltarr(nx,ny,nt)
spike=fltarr(nx,ny,nt)
FOR i=0,nt-1 DO BEGIN

   dum=cube[*,*,i]
   atrous,dum,decomp=a,n_scales=1
   despike[*,*,i]=a[*,*,0]
   spike[*,*,i]=a[*,*,1]
ENDFOR

FOR i=0,nt-1 DO BEGIN
   z=spike[*,*,i]
   m=moment(z[where(dum ne 0)])
   in=where(z gt 3.*sqrt(m[1]))
   ;z2=gauss_smooth(z,2,/edge_truncate)
   ;in=where(z ne 0 and r gt 233)
   ;res=moment(z2[in])
   ;in2=where(abs(z[in]-res[0]) gt 5.*sqrt(res[1]))
   ;z[in[in2]]=0.
   ;z3=gauss_smooth(z,3,/edge_truncate) 
   ;z[in[in2]]=z3[in[in2]]
   IF i lt 2 THEN t1=0 ELSE t1=i-2
   IF i gt nt-3 THEN t2=nt-1 ELSE t2=i+2
   mn=median(spike[*,*,t1:t2],dim=3)
   mn=median(mn,5)
   z[in]=mn(in)
   despike[*,*,i]+=z  


ENDFOR

stop



END