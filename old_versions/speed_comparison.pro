pro speed_comparison, zoom=zoom, noscl=noscl, scl=scl, force_zero=force_zero

;wave angle sav-file written by wave_tracking.pro
restore, 'wave_angle_small_20120327_3_3.5mHz.sav'
;speeds asav-file written by wave_tracking.pro if /wr_speeds-keyword was set
restore, 'speeds_testing_small_20120327_3_3.5mHz.sav'

data = reform(wave_angle)
pbla = [0,0,0,0] 
rbla = [0,0,0,0] 

if not keyword_set(zoom) then zoom = 1.
dim = size(data)
sx = dim[1]
sy = dim[2]
nsx = sx*zoom
nsy = sy*zoom
if (dim[0] ne 2) then begin
   print, 'Sorry - data must be a 2D-array'
   GOTO, CHEERIO
endif
if (nsx lt 350) then begin
   print,'****************************************************'
   print,'* Please use at least 350 pixels for a proper view *'
   print,'****************************************************'
   print,''
endif

print, '*********************************************'
print, '* Press right mouse button in image to quit *'
print, '*********************************************'

;=== initial definitions ===
!MOUSE.BUTTON = 0
x = 0. & xstr = ' '
y = 0. & ystr = ' '
lastx = -1.
lasty = -1.
value = ' '
new_data = congrid(data, nsx, nsy, cubic=0)
if (n_elements(scl) eq 2) then new_data = bytscl(new_data,min=scl[0],max=scl[1])

window, xs=nsx+100., ys=nsy+100., /free
new_w = !D.WINDOW
window, xs=900, ys=600, /free
newer_w = !D.WINDOW
wset, new_w

loadct, 4, /silent

if (keyword_set(noscl) or (n_elements(scl) eq 2)) then begin
   tv, new_data, 50., 50.
endif else begin
   tvscl, new_data, 50., 50.
endelse

while (!MOUSE.BUTTON ne 4) do begin
 cursor, x, y, 2, /dev
 xyouts, nsx-40.,  nsy+82., 'x: '+xstr,      charsize=2., /dev, col=0. ; erase old print
 xyouts, nsx-39.,  nsy+63., 'y: '+ystr,      charsize=2., /dev, col=0. ; ===============
 wset, newer_w
 xyouts, 80, 580, 'New prograde speed: '+strcompress(string(pbla[0]*1000.,format='(F12.2)'),/rem)+' '+string(177b)+' '+strcompress(string(pbla[1]*1000.,format='(F8.2)'),/rem)+' km/s', charsize=1.5, color=0, /dev
 xyouts, 80,  280, 'New retrograde speed: '+strcompress(string(rbla[0]*1000.,format='(F12.2)'),/rem)+' '+string(177b)+' '+strcompress(string(rbla[1]*1000.,format='(F8.2)'),/rem)+' km/s', charsize=1.5, color=0, /dev
 xyouts, 520, 580, 'Old prograde speed: '+strcompress(string(pbla[2]*1000.,format='(F12.2)'),/rem)+' '+string(177b)+' '+strcompress(string(pbla[3]*1000.,format='(F8.2)'),/rem)+' km/s', charsize=1.5, color=0, /dev
 xyouts, 520, 280, 'Old retrograde speed: '+strcompress(string(rbla[2]*1000.,format='(F12.2)'),/rem)+' '+string(177b)+' '+strcompress(string(rbla[3]*1000.,format='(F8.2)'),/rem)+' km/s', charsize=1.5, color=0, /dev
 wset, new_w  
 x = fix((x-50.)/zoom)
 y = fix((y-50.)/zoom)
 device, decomposed=0

 if (keyword_set(noscl) or (n_elements(scl) eq 2)) then begin
    tv, new_data, 50., 50.
  endif else begin
    tvscl, new_data, 50., 50.
 endelse

 xstr = strtrim(strcompress(string(x)),2) & ystr = strtrim(strcompress(string(y)),2)
 if ((0 gt x) or (x gt sx-1.) or (0 gt y) or (y gt sy-1.)) then begin 
    xstr = 'no data'
    ystr = 'no data'
 endif else begin
    loadct, 39, /silent
    !p.multi = [0,2,2]
    wset, newer_w
    if keyword_set(force_zero) then begin
     pbla=compute_speed_comparison(reform(p_vel[x,y,*,*]),3.23924,30,debug=1,/force_zero)
     rbla=compute_speed_comparison(reform(r_vel[x,y,*,*]),3.23924,30,debug=1,/ret,/force_zero)
    endif else begin
     pbla=compute_speed_comparison(reform(p_vel[x,y,*,*]),3.23924,30,debug=1)
     rbla=compute_speed_comparison(reform(r_vel[x,y,*,*]),3.23924,30,debug=1,/ret)
    endelse
    xyouts,80, 580, 'New prograde speed: '+strcompress(string(pbla[0]*1000.,format='(F12.2)'),/rem)+' '+string(177b)+' '+strcompress(string(pbla[1]*1000.,format='(F8.2)'),/rem)+' km/s', charsize=1.5, color=254, /dev
    xyouts, 80,  280, 'New retrograde speed: '+strcompress(string(rbla[0]*1000.,format='(F12.2)'),/rem)+' '+string(177b)+' '+strcompress(string(rbla[1]*1000.,format='(F8.2)'),/rem)+' km/s', charsize=1.5, color=254, /dev
    xyouts, 520, 580, 'Old prograde speed: '+strcompress(string(pbla[2]*1000.,format='(F12.2)'),/rem)+' '+string(177b)+' '+strcompress(string(pbla[3]*1000.,format='(F8.2)'),/rem)+' km/s', charsize=1.5, color=254, /dev
    xyouts, 520, 280, 'Old retrograde speed: '+strcompress(string(rbla[2]*1000.,format='(F12.2)'),/rem)+' '+string(177b)+' '+strcompress(string(rbla[3]*1000.,format='(F8.2)'),/rem)+' km/s', charsize=1.5, color=254, /dev
    wset, new_w
    loadct, 4, /silent
 endelse

 xyouts, nsx-40., nsy+82., 'x: '+xstr,  charsize=2., color=254, /dev
 xyouts, nsx-39., nsy+63., 'y: '+ystr,  charsize=2., color=254, /dev

endwhile

if (!D.WINDOW eq new_w) then wdelete, new_w & wdelete, newer_w

CHEERIO:
stop
end
