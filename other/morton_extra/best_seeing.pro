pro best_seeing

fname='analysis/comp/seeing.txt'
fmt='L,F,F' 
readcol,fname,date,time,fwhm,format=fmt,/quick,count=count

h=0


WHILE h lt count-1 do begin
  
    sdat=date[h]
    in=where(date eq sdat,num)
    fwhm_date=fwhm[in]
    in2=where(fwhm_date lt 50) 
    IF n_elements(avsee) EQ 0 THEN avsee=[(moment(fwhm_date[in2]))[0:1],median(fwhm_date[in2])] $    
    ELSE avsee=[[temporary(avsee)],[[(moment(fwhm_date[in2]))[0:1],median(fwhm_date[in2])]]]

    IF n_elements(date2) EQ 0 THEN date2=sdat $    
    ELSE date2=[temporary(date2),sdat]
    
    h=h+num
    
ENDWHILE

pdf = HISTOGRAM(avsee[0,*], LOCATIONS=xbin,binsize=0.1)
p=plot(xbin,pdf, title='Mean values',layout=[2,1,1])

pdf = HISTOGRAM(avsee[2,*], LOCATIONS=xbin,binsize=0.1)
p2=plot(xbin,pdf, title='Median values',layout=[2,1,2],/current)

mn=moment(avsee[0,*])
print,'Moment of average seeing', mn

in=where(avsee[0,*] lt mn[0]-sqrt(mn[1]))
print,'Days of best seeing',date2[in]

;date2s=strtrim(date2,2)
;year=strmid(date2s,0,2)
;month=strmid(date2s,2,2)
;day=strmid(date2s,4,2)
;juldate=julday(month,day,year)

;dummy = LABEL_DATE(DATE_FORMAT=['%D-%M','%Y'])
;dplot=plot(juldate[in],avsee[0,in],xtickunits = ['Time', 'Time'], $
;   xtickformat='LABEL_DATE', XSTYLE=1 )
stop

END