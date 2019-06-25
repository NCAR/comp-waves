pro progress_comp, cur, max, msg=msg
  if n_elements(msg) eq 0 then begin
   msg=""
   writeu, -1, string(format='(%"\R",A,"",I3,"% done.")', $
               msg, round(float(cur)/max*100))
  endif else begin 
   writeu, -1, string(format='(%"\R",A,": ",I3,"% done.")', $
               msg, round(float(cur)/max*100))
  endelse
  if cur eq max then print, ''
end
