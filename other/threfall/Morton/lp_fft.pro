pro lp_fft,file ,t0=t0,t1=t1
;+
;   lp_fft,file ,t0=t0,t1=t1
;
;            calculates fft of cube and writes abs(power) to file_fft.fcube
;            works in core so memory demanding
;-
if(n_params() lt 1) then begin
  message,'calling sequence: lp_fft,file ,t0=t0,t1=t1',/info
  return
endif
dum=file_search(file,count=count)
if(count eq 0) then begin
  message,'file not found',/info
  return
endif
ip=strpos(file,'.icube')
if(ip lt 0) then ip=strpos(file,'.bcube')
if(ip lt 0) then begin
  print,'file has to be type .bcube or .icube'
  return
endif

file_base=strmid(file,0,ip)
file_fft=file_base+'_fft_xfy.fcube'

lp_header,file,nx=nx,ny=ny,nt=nt
if(n_elements(t0) eq 0) then t0=0
if(n_elements(t1) eq 0) then t1=nt-1
extraheader='t0='+strtrim(string(t0),2)+', t1='+strtrim(string(t1),2)
im=lp_read(file)
i=0
pow=abs(fft(reform(im[*,i,t0:t1]),-1,dimension=2))^2
lp_put,pow[*,0:(t1-t0+1)/2],file_fft,0,nt=ny,/keep_open,extraheader=extraheader
print,''
for i=1,ny-1 do begin
  pow=abs(fft(reform(im[*,i,t0:t1]),-1,dimension=2))^2
  lp_put,pow[*,0:(t1-t0+1)/2],file_fft,i,/keep_open
  print,string(13b)+'% finished: ',float(i)*100./(ny-1),format='(a,f4.0,$)'
endfor
print,''
print,'transposing'
pow=0
im=0
im=lp_read(file_fft)
im=temporary(transpose(im,[0,2,1]))
file_fft=file_base+'_fft_xyf.fcube'
lp_write,im,file_fft,extraheader=extraheader

end
