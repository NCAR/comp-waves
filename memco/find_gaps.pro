;+
;  Name: find_gaps
;
;  Description:
;    Procedure to identify gaps in a time series. Pass this routine a one dimensional array of data points
;    where missing data points are set to 0. 
;
;  Input:
;    data - one dimensional array of data 
;
;  Output:
;    gap_start - an array of the starting indices of the gaps (of size number of gaps)
;    gap_end - an array of the ending indices of the gaps (of size number of gaps)
;    ngap - number of gaps identified
;
;  Keyword Parameters: none
;
;  Author: Tomczyk
;
;  Examples:
;    find_gaps, data, gap_start, gap_end, ngap
;-
pro find_gaps, data, gap_start, gap_end, ngap
 
;  procedure to find gaps in data strings, identified by zero data
 
debug='yes'		;debug mode ('yes' or 'no')

gap_start=intarr(5000)    ;dimension to largest number of possible gaps
gap_end=intarr(5000)
num=n_elements(data)
 
ngap=0
time=where( data ne 0., count)
if count gt 0 then begin
 
for i=1,count-1 do begin
  if time(i)-time(i-1) gt 1 then begin
    if debug eq 'yes' then begin
      print, time(i-1), time(i)
      print, 'gap at i=',i,' at times',time(i-1)+1,' to',time(i)-1,' inclusive'
    endif
    gap_start(ngap)=time(i-1)+1
    gap_end(ngap)=time(i)-1
    ngap=ngap+1
  endif
endfor

if ngap gt 0 then begin
  gap_start=gap_start(0:ngap-1)
  gap_end=gap_end(0:ngap-1)
endif

endif

end
