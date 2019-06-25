pro lp_put, image, filename, indx, nt=nt, extraheader=extraheader, $
 swap_endian=swap_endian, keep_open=keep_open, close=close,verbose=verbose
;
; writes La Palma data as IDL assoc file
; first 512 bytes are used to store header info
; this routine writes one image at a time
; header is written if nt keyword is given
; minimum required header info for lp_read is then written here, additional
; header info can be added in <extraheader>
; /keep_open does not close the file and assumes it is open if indx > 0
;
common clp_put,luw,indx_mx

 IF n_params() LT 3 THEN BEGIN
   message, /info, 'lp_put, image, filename, indx [, extraheader=extraheader, /swap_endian]'
   return
 ENDIF

 if KEYWORD_SET(extraheader) eq 0 then extraheader=''
 if KEYWORD_SET(swap_endian) eq 0 then swap_endian=0

 sZ = size(image)
 dims = sZ[0]
 if dims ne 2 then begin
   message, /info, ' only 2D image is supported'
   return
 endif
 nx = sZ[1]
 ny = sZ[2]
 datatype = sZ[dims+1]

 files=file_search(filename,count=count)
 if count EQ 0 THEN BEGIN
   openw, luw, filename, /get_lun, swap_endian=swap_endian
 endif else begin
   if(not keyword_set(keep_open)) or (indx eq 0) then begin
     openu, luw, filename, /get_lun, swap_endian=swap_endian
   endif
 endelse
 if n_elements(nt) NE 0 THEN BEGIN
   indx_mx=nt-1
   bheader = make_lp_header(image)
   ipos=strpos(bheader,'dims=')
   strput,bheader,'3',ipos+5
   ipos=strpos(bheader,', endian=')
   text1=strmid(bheader,0,ipos)
   text2=strmid(bheader,ipos,strlen(bheader)-ipos)
   if swap_endian NE 0 THEN BEGIN
     endian = strmid(text2,strlen(text2)-1,1)
     if endian EQ 'b' THEN endian='l' ELSE endian='b'
     text2=strmid(text2,0,strlen(text2)-1)+endian
   endif
   bheader=text1+', nt='+strtrim(string(nt),2)+text2
   if(keyword_set(verbose)) then message, /info, bheader
   header = extraheader +' : ' + bheader
   rec = assoc(luw, bytarr(512))
   rec[0] = byte(header)
 endif
 case datatype of
     1: begin
         rec = assoc(luw, bytarr(nx,ny), 512)
         rec[indx] = image
     end
     2: begin
         rec = assoc(luw, intarr(nx,ny), 512)
         rec[indx] = image
     end
     3: begin
         rec = assoc(luw, lonarr(nx,ny), 512)
         rec[indx] = image
     end
     4: begin
         rec = assoc(luw, fltarr(nx,ny), 512)
         rec[indx] = image
     end
     else : begin
         message, ' datatype not supported'
         print, ' datatype = ', datatype
     end
 endcase
 if(not keyword_set(keep_open)) then begin
   free_lun,luw
 endif else if (indx eq indx_mx) then free_lun, luw

end
