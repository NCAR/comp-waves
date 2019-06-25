; datadir = '/mn/astro-cluster/u1/lapalma/soup_02Jun2003/mfbd/series/'
; cubename = 'ni6768_RCP_strt.icube'
; cubename = 'ni6768_magneto.icube'
; lp_ximovie, datadir+cubename, mag=.4
; datadir='/mn/astro-cluster/u1/lapalma/gband_02Jun2003/mfbd/series/'
; cubename = 'gband_02Jun2003_strt.icube'
; lp_ximovie, datadir+cubename, mag=.4

pro lp_ximovie, cubename, magnification=magnification,_extra=xikeywords
;+
; start ximovie.pro on La Palma data cube
; reads info from header
;-
 if n_params() eq 0 then begin
     message, /info, 'lp_ximovie, cubename [, ximovie_keywords]'
     return
 endif

 hoffset = 512   ; header is 512 bytes
 lp_header, cubename, datatype=datatype, nx=nx, ny=ny, nt=nt, endian=endian_file
 if(n_elements(endian_file) eq 0) then endian_file='l'
 if ((byte(1L, 0, 1))[0] eq 1) then endian = 'l' else endian='b'
 if(datatype gt 1) and (endian ne endian_file) then swap_endian=1 else swap_endian=0
 if(swap_endian eq 1) then begin
   message,/info,'swapping from endian='+endian_file+' to endian='+endian
 endif
; make sure window fits on screen
 if(n_elements(magnification) eq 0) then begin
   magnification=1.0
   device,get_screen_size=screen_size
   npx_decoration=[281,95]   ; number of pixels for sliders etc
   if(screen_size[0] gt 2000) then screen_size[0]=screen_size[0]/2.
   max_window=screen_size-npx_decoration
   if(nx gt max_window[0]) then magnification=float(max_window[0])/nx
   if(ny gt max_window[1]) then magnification=magnification < float(max_window[1])/ny
 endif
 cubesize = float(nx) * float(ny) * float(nt) * float(datatype)
 case datatype of
     1: typestring='bytarr'
     2: typestring='intarr'
     4: typestring='fltarr'
     else: typestring='unsupported type: '+strtrim(datatype,2)
 endcase
 si = ' cube is ('+strtrim(nx,2)+', '+strtrim(ny,2)+', '+strtrim(nt,2)+'), '+typestring
 si = si + ' = '+strtrim(cubesize/1.e6,2)+' Mb'
 print, si
 print, '  starting ximovie...'
 case datatype of 
     1 : ximovie, cubename, nx, ny, offset=hoffset, magnification=magnification,_extra=xikeywords
     2 : ximovie, cubename, nx, ny, offset=hoffset, magnification=magnification,_extra=xikeywords, /int, swap_endian=swap_endian
     4 : ximovie, cubename, nx, ny, offset=hoffset, magnification=magnification,_extra=xikeywords, /float, swap_endian=swap_endian
     else : print, 'datatype not supported by ximovie, datatype='+strtrim(datatype,2)
 endcase

 return
end
