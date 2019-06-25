;+
;  Name: fill_gaps
;
;  Description:
;    Procedure to fill gaps in a data series using the maximum entropy method. (Fahlman and Ulrych, 9182).
;    This routine calls the routine memcof from Numerical Recipes. First call the find_gaps routine. Then  
;    pass this routine a one dimensional array of data points where missing data points are set to 0 adata
;    with gap_start, gap_end, and ngap output from find_gaps. 
;
;  Notes: If distance between gaps is small and coef_num is larger, then to keep coef_num the same
;         the current version will estimates the LP coefficients will already filled data
;
;  Input:
;    data - one dimensional array of data
;    coef_num - number of Linear Prediction coefficients to estimate
;
;  Output:
;    data - the gap filled data array is returned
;
;  Keyword Parameters: none
;
;  Author: Tomczyk
;  Edits: RJ Morton APR2019 - removed loops. Corrected implementation. Previous version used entire time-series to 
;                     to estimate LP coefficients. Should only use non-gapped series before and after
;                     
;
;  Example:
;    find_gaps, data, gap_start, gap_end, ngap
;    fill_gaps, data, gap_start, gap_end, ngap
;-
pro fill_gaps, data, gap_start, gap_end, ngap,coef_num=coef_num,verbose=verbose,aic_sel=aic_sel,$
                types=types

;IF n_elements(npoles) EQ 0 THEN npoles=35    ;the number of poles can be adjusted
IF n_elements(coef_num) EQ 0 THEN coef_num=35 ;Number of coefficients - M

debug="no"   ;debug mode ('yes' or 'no')

max_gap=40      ;maximum gap length to fill (can be adjusted)

;  find maximum entropy coefficients
data_size=n_elements(data)
no_zero=where(data ne 0., count)
first=min(no_zero)	    ;first real data point
last=max(no_zero)	      ;last real data point


gap_size=gap_end-gap_start+1

;Find longest stretch of time-series to calculate LP coefficients
;Inherently assume stationarity of series
IF ngap gt 1 THEN BEGIN
 data_cont=gap_start[1:-1]-gap_end[0:-2]
 data_cont=[gap_start[0],data_cont,data_size-gap_end[-1]]
 
ENDIF ELSE BEGIN
  data_cont=[gap_start[0],data_size-gap_end[0]]
ENDELSE

mx=max(data_cont,loc)
if mx lt coef_num THEN coef_num=mx ELSE $
                                     coef_num=coef_num
IF loc eq 0 THEN BEGIN 
    t1=0
    t2=gap_start[0]-1
ENDIF ELSE BEGIN
    IF ngap eq loc THEN BEGIN
       t1=gap_end[-1]+1
       t2=data_size-1
    ENDIF ELSE BEGIN
        t1=gap_end[loc-1]+1
        t2=gap_start[loc]-1        
    ENDELSE
ENDELSE


;Calculates Linear prediction coefficients


; IF keyword_set(aic_sel) THEN BEGIN
;      cont_len=t2-t1 ;number of data points in longest stretch
;      xm=fltarr(40)
;      min_val=10
;      max_val=50<cont_len

;      FOR i=min_val,max_val DO BEGIN
        if memcof(data(t1:t2), coef_num, xms, coef) ne 0 then begin
            print, 'error in memcof'
        endif
;        IF keyword_set(verbose) THEN print,'Mean Square Discrepancy',xms,i
;        xm[i-min_val-1]=xms
;      ENDFOR
;ENDIF

;Fixes bad guesses for coefficients - applies a stability condition
fixrts, coef, coef_num

types=fltarr(ngap)
for i=0,ngap-1 do begin
    gst=gap_start[i]
    gend=gap_end[i]
    input=dblarr(coef_num)
    
    type=0  ;type of filling (0=no, 1=forward, 2=backward, 3=both) 
   
    ;  fill gaps 
    if gap_size[i] le max_gap then begin	;only fill gaps le max_gap samples
         
        ;  predict forward
        forward_pre=dindgen(coef_num)-coef_num + gst
        good=where( data(forward_pre) ne 0., ngood )
        
        ;verify that predictors are real data
        ;Condition checks no gaps exist before and whether enough data points
        ;available before hand
        if ngood eq coef_num and gst-coef_num ge 0 then begin
            input=data[forward_pre]
           ;for j=0,coef_num-1 do input(j)=data(gst-coef_num+j)
            future=predic(input,coef,gap_size[i])
            type=type+1
        endif
        if debug eq 'yes' then if ngood ne coef_num then print, ngood
        
        ;  predict backward     
        back_pre=dindgen(coef_num) + gend+1
        good=where( data(back_pre) ne 0., ngood )
        
        ;  verify that predictors are real data
        ;Condition checks no gaps exist after and whether enough data points
        ;available
        if debug eq 'yes' then if ngood ne coef_num then print, ngood,coef_num,gend
        
        if ngood eq coef_num and gend+coef_num lt data_size then begin   
            input=data[back_pre]
            ;for j=0,coef_num-1 do input(j)=data(gend+coef_num-j)
            future2=predic(reverse(input),coef,gap_size[i])

            future2=reverse(future2)
            type=type+2
        endif     
       
        if debug eq 'yes' then print,type

        case type of
            0: ;print,'gap could not be filled'
            1: data(gst:gend)=future
            2: data(gst:gend)=future2 
            3: data(gst:gend)=(future+future2)/2. ;  average forward and backward predictions
        endcase     
    endif
    types[i]=type
endfor
 
end
