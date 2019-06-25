pro blue_red, low, top, center=center, width=width, slope=slope, mix=mix
;+
; NAME:
;
;   red_blue
;-
N_colors=!D.n_colors<256

IF not keyword_set(low)    then low=0

IF not keyword_set(top)    then high=N_colors-1 $
  			   else high=N_colors-1-top

If not keyword_set(slope)  then slope=1. ELSE slope=float(slope)>0.01

If not keyword_set(center) then center = byte((low+high)/2.)

IF not keyword_set(mix)    then mix=0

IF not keyword_set(width)  then width=0

;
; load the definitions of the colors
;  
tvlct, r, g, b, /get

r=bytarr(N_colors)
g=bytarr(N_colors)
b=bytarr(N_colors)

inc=bytscl( (findgen(center-low+1)+1) ^ (1./float(slope)) )

dec=reverse(bytscl( (findgen(high-center+1)+1) ^ (1./float(slope)) ))

if keyword_set(inverse) then $
     begin
     r(low:center)=255  &  r(center:high)=dec
     g(low:center)=inc  &  g(center:high)=dec
     b(low:center)=inc  &  b(center:high)=255
     end $
   else $
     begin
     r(low:center)=inc  &  r(center:high)=255
     g(low:center)=inc  &  g(center:high)=dec
     b(low:center)=255  &  b(center:high)=dec
     end

If keyword_set(verbose) then $
    begin
    p=!p.multi & !p.multi=[0,1,3]
    plot,r,xs=3,ys=3,chars=2&plot,g,xs=3,ys=3,chars=2&plot,b,xs=3,ys=3,chars=2
    !p.multi=p
    end

IF high LT N_colors-1 THEN $
   BEGIN
   r(high+1:N_colors-1)=255
   b(high+1:N_colors-1)=255
   g(high+1:N_colors-1)=255
   END

;
; store the definitions of the colors
;  
tvlct, r, g, b	; save new color table

; !P.multi=[0,1,3]
; plot,r,xs=1,chars=2
; plot,g,xs=1,chars=2
; plot,b,xs=1,chars=2
; !P.multi=0

end 
