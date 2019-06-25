pro display_kimage,image,mn,mx,xt=xt,yt=yt,tit=tit,units=units,table=table,$
 output=output,windnum=windnum,wherebar=wherebar,scale=scale,offset=offset,$
 phase_speed=phase_speed

;  Subroutine to display k-omega diagram.

;  table = color table number
;  units = units for color bar
;  mn,mx = min and max for image scaling and colorbar limits
;  output = output device 'p' or 's'
;  wherebar = location of color bar 'top' or 'right'

ans=' '

s=size(image)
nx=s(1)
ny=s(2)

if n_elements(offset) eq 0. then offset=0.

img=bytscl(image,mn,mx)

;  set margins for plot (in pixels)

if n_elements(yt) eq 0. then left=25. else left=80.
if wherebar eq 'right' then begin
  bottom=80.
  right=200.
  top=50.
endif
if wherebar eq 'top' then begin
  bottom=80.
  right=25.
  top=120.
endif
if wherebar eq 'bottom' then begin
  bottom=150.
  right=25.
  top=50.
endif


;  compute size of window in pixels

xw=float(nx)+left+right
yw=float(ny)+top+bottom

;  open window if screen display

if output ne 'p' then begin
  window,windnum,xs=xw,ys=yw
  device,decomposed=0
endif

;  load color table, insure white border

loadct,table, /silent
v=intarr(256,3)
tvlct,v,/get
v(255,*)=255
tvlct,v
img=img<254

xmin=offset
xmax=scale*float(nx-1)+offset
ymin=offset
ymax=scale*float(ny-1)+offset

;  display image

tv,img,left,bottom,xsize=nx,ysize=ny

tny=1./(2.*28.7)        ;temporal nyquist (29 s sampling)
sny=1./(2.*3.22)       ;spatial nyquist (3.22 Mm/pixel)

plot,[-tny,tny],[0.,sny],/nodata,chars=1.5,title=tit, $
 position=[left,bottom,left+nx,bottom+ny],xsty=1,/device, $
 xtitle=xt,ytitle=yt,/noerase,ystyle=1,ticklen=0.01

xyouts,.2,.5,'Upward',/norm,chars=2
xyouts,.55,.5,'Downward',/norm,chars=2

if n_elements(phase_speed) eq 0. then phase_speed=0.6  ;constant phase velocity (Mm/s) (0.6=default)
x=tny*findgen(100)/101.
y=x/phase_speed
oplot,x,y,linesty=1
oplot,-x,y,linesty=1


;  draw vertical color bar on right

if wherebar eq 'right' then begin
  xc=16          ;x size
  yc=fix(ny*0.8)

  xpos=left+nx+right/2.
  ypos=bottom+0.1*ny

  col=254.*findgen(yc)/float(yc-1)  ;244 to omit white at top of color table
  col=rebin(col,yc,xc)
  col=transpose(col)
  tv,col,xpos,ypos

  plot,[0,1],[mn,mx],/nodata,chars=1.5, $
   position=[xpos,ypos,xpos+xc,ypos+yc],/device,xsty=1,xticks=1, $
   ytitle=units, /noerase, ystyle=1, ticklen=0.04,xchars=1.e-4
endif

;  draw horizontal color bar on top

if wherebar eq 'top' then begin
  yc=16          ;x size
  xc=fix(nx*0.8)

  xpos=fix(left+0.1*nx)
  ypos=fix(bottom+ny+0.6*top)

  col=254.*findgen(xc)/float(xc-1)
  col=rebin(col,xc,yc)
  tv,col,xpos,ypos

  plot,[mn,mx],[0,1],/nodata,chars=1.5, $
   position=[xpos,ypos,xpos+xc,ypos+yc],/device,xstyle=1,yticks=1, $
   xtitle=units, /noerase, ystyle=1, ticklen=0.04,ychars=1.e-4
endif

;  draw horizontal color bar on bottom

if wherebar eq 'bottom' then begin
  yc=16          ;x size
  xc=fix(nx*0.8)

  xpos=fix(left+0.1*nx)
  ypos=fix(0.4*bottom)

  col=254.*findgen(xc)/float(xc-1)
  col=rebin(col,xc,yc)
  tv,col,xpos,ypos

  plot,[mn,mx],[0,1],/nodata,chars=1.5, $
   position=[xpos,ypos,xpos+xc,ypos+yc],/device,xstyle=1,yticks=1, $
   xtitle=units, /noerase, ystyle=1, ticklen=0.04,ychars=1.e-4
endif

end
