FUNCTION cutter, data, x1, x2, y1, y2
; program to calculate the pixel values in a cut across data (x1,y1)-(x2,y2)
; data is a 2d pixel window
; cut begins at x1,y1 and ends at x2, y2
; NB 2d functionality still questionable - 1d cuts have been tested

nx=abs(x2-x1)
ny=abs(y2-y1)

;figure out m and c for line y=mx+c
m=double(y2-y1)/double(x2-x1)
c=y1-m*x1

;if x or y is constant, then we simply count each pixel along the line
 IF ((ny gt nx) and (nx eq 0)) THEN BEGIN
  intcoords=dblarr(ny+1,3)
  FOR i=0,ny-1 DO BEGIN
   intcoords(i,0)=x1
   IF (y2 gt y1) THEN intcoords(i,1)=y1+i
   IF (y1 gt y2) THEN intcoords(i,1)=y2+i
   IF (y2 gt y1) THEN intcoords(i,2)=data(x1,y1+i)
   IF (y1 gt y2) THEN intcoords(i,2)=data(x1,y2+i)
  ENDFOR
 ENDIF
 IF ((nx gt ny) and (ny eq 0)) THEN BEGIN
  intcoords=dblarr(nx+1,3)
  FOR i=0,nx-1 DO BEGIN
   IF (x2 gt x1) THEN intcoords(i,1)=x1+i
   IF (x1 gt x2) THEN intcoords(i,1)=x2+i
   ;intcoords(i,0)=x1+i
   intcoords(i,1)=y1
   IF (x2 gt x1) THEN intcoords(i,2)=data(x1+i,y1)
   IF (x1 gt x2) THEN intcoords(i,2)=data(x2+i,y1)
  ENDFOR
 ENDIF
;if we move in more than 1d, count pixels in longer axis and use those values.
 IF ((ny gt nx) and (nx ne 0)) THEN BEGIN
  intcoords=dblarr(ny+1,3)
  IF (y2 gt y1) THEN  xvals=ROUND(((y1+findgen(ny+1))-c)/m)
  IF (y1 gt y2) THEN  xvals=ROUND(((y1-findgen(ny+1))-c)/m)
  FOR i=0,ny-1 DO BEGIN
   intcoords(i,0)=xvals(i)
   IF (y2 gt y1) THEN BEGIN  ; if y2 biggest, start at y1 and work to y2
   intcoords(i,1)=y1+i
   intcoords(i,2)=data(xvals(i),y1+i)
   ENDIF ELSE BEGIN
   ;IF (y1 gt y2) THEN BEGIN	;if y1 biggest, start at y1 and go down!
   intcoords(i,1)=y1-i
   intcoords(i,2)=data(xvals(i),y1-i)
   ENDELSE
  ENDFOR
 ENDIF
 IF ((nx ge ny) and (ny ne 0)) THEN BEGIN 
  intcoords=dblarr(nx+1,3)
  IF (x2 gt x1) THEN yvals=ROUND(m*(x1+findgen(nx+1))+c)
  IF (x1 gt x2) THEN yvals=ROUND(m*(x1-findgen(nx+1))+c)
  FOR i=0,nx-1 DO BEGIN
   intcoords(i,1)=yvals(i)
   IF (x2 gt x1) THEN BEGIN
   intcoords(i,0)=x1+i
   intcoords(i,2)=data(x1+i,yvals(i))
   ENDIF ELSE BEGIN
   intcoords(i,0)=x1-i
   intcoords(i,2)=data(x1-i,yvals(i))
   ENDELSE
  ENDFOR 
 ENDIF


;wait, 0.5
radialf=dblarr(n_elements(intcoords(*,0))-1,2)

;store the resulting data in terms using a radial distance and f(r)
FOR i=0,n_elements(intcoords(*,0))-2 DO BEGIN
 ;oplot, [icoords(i,0),icoords(i,0)], [icoords(i,1),icoords(i,1)], psym=5, col=180;, symsize=2

radialf(i,0)=sqrt((intcoords(i+1,0)-x1)*(intcoords(i+1,0)-x1)+(intcoords(i+1,1)-x1)*(intcoords(i+1,1)-x1))
 radialf(i,1)=intcoords(i+1,2)
ENDFOR
;plot, radialf(*,1)
;STOP
;return the radial data
RETURN, radialf

END
