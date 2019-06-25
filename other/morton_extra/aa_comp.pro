;PURPOSE: - Calculate angular-ly averaged velocity power spectra for CoMP data
;
;
;METHODOLOGY:
;           The aim is to fit the power spectra with a function P=10^af^-b. This will be two
;           parts due to oscillatory bump in power spectra.
;           
;           Using FFT to find periodogram. Following Vaughan, A&A, 431,391 (2005) the routine
;           finds FFT power and converts to log space to make errors homoskedastic.
;           In this case, power in each frequency bin should be close to normal distribution (assuming
;           the sample of the corona we choose all has the same power spectra).
;
;           For each frequency bin, the mean log power and standard error is calculated through 
;           the bootstrap (Enfron) and a Kolmogorov-Smirnoff test is run to check for normality
;           of the distribution of the mean log power. This is done in power_fit_cube. Details of 
;           each frequency distribution are also used.
;
;            
; 
;
;
;INPUTS:- date  - date of aligned cube_ivw
;
;OPTIONAL INPUTS: - apod - degree of apodisation (<1.) - default is 0.1
;                 - extra - string relating to additional characters in cube_ivw name
;                 - /doplot - Outputs diagnostic plots
;                 - phi_ang - select the angle to average over
;                 - rmin - define a minimum radius for averaging
;                 - rmax - define a maximum radius for averaging
;
; EXAMPLE - aa_comp,'20120327',phi_ang=5,rmax=1080,rmin=1050   
;
;
;OUTPUTS: - 
;
;TO DO: - 
;         Problem with 'duff' time-series - current removal of time-series with number of zeros is problematic
;         Might be better to exclude series on rms power values
;


;----------------------------------------------------------------------------
;Calculate r,phi values for image using coord_cart_heilo
function frad_coord,index,image,phi=phi

sz=size(image)                                          
nx=sz(1)
ny=sz(2)      
x=findgen(nx)
x=rebin(x,nx,ny)
y=findgen(ny)   
y=rebin(y,nx,ny)
y=transpose(y)   

coord_cart_helio,index,1.0,x,y,xout,yout,l,b

r=sqrt(xout^2+yout^2)


x_in=(where(xout[*,0]*shift(xout[*,0],1) lt -1))[1]
y_in=(where(xout[*,0]*shift(xout[*,0],1) lt -1))[1]

phi=fltarr(nx,ny)
phi=90-atan(yout/xout)*180./!pi

phi[0:x_in-1,0:y_in-1]=270.-atan(yout[0:x_in-1,0:y_in-1]/xout[0:x_in-1,0:y_in-1])*180./!pi 
phi[0:x_in-1,y_in:ny-1]=270.-atan(yout[0:x_in-1,y_in:ny-1]/xout[0:x_in-1,y_in:ny-1])*180./!pi 

return,r

END

;----------------------------------------------------------------------------




;----------------------------------------------------------------------------
PRO AA_comp,date,index=index,$           
                extra=extra,doplot=doplot,$
                phi_ang=phi_ang,rmin=rmin,rmax=rmax


; paths where the CoMP data is located and
; where the .sav-files should be written to
inpath  = 'analysis/CoMP/wave_tracking_output/'+date+'/'
outpath = 'analysis/CoMP/wave_tracking_output/'+date+'/'



IF NOT keyword_set(extra) THEN restore,inpath+'cube_ivw_'+date+'.sav',/verbose $
ELSE restore,inpath+'cube_ivw_'+date+'_'+extra+'.sav',/verbose

data=cube_w
sz=size(data)
nx=sz(1) & ny=sz(2) & nt=sz(3)
xscale=index.cdelt1 ;pixel size in arcsec


; Upper and lower limit for the radius in
; !be aware that the maximum offset of the cross-correlation
; (if performed) is added to this value!
IF n_elements(rmin) EQ 0 THEN lower_r = index.lower_r*xscale ELSE lower_r=rmin
IF n_elements(rmax) EQ 0 THEN upper_r = index.upper_r*xscale ELSE upper_r=rmax
maxoffset=index.maxoffset*xscale

; Create arrays of r & phi values
r=frad_coord(index,data[*,*,0],phi=phi)

; create mask & mask out pixels not to invert
mask=index.mask
mask[where(r lt lower_r)]=0
mask[where(r gt upper_r)]=0
data_mask=data
for i=0,nt-1 do data_mask[*,*,i]=data[*,*,i]*mask


;Phi averaging
IF n_elements(phi_ang) EQ 0 THEN phi_ang=45
no_phi=360/phi_ang

;data store
temp={date:index.date_d$obs,phi:0.0,mn:0.0,med:0.0,sd:0.0}


print,'Starting angle averaging'
FOR i=1,no_phi-1 DO BEGIN
     phold=!null
     rval=!null

     ;Determine pixels that lie between two angles          
     wedge=where(phi gt phi_ang*(i-1) and phi lt phi_ang*i and mask gt 0)
     xc=wedge mod nx
     yc=wedge/nx
     
     
     FOR j=0,n_elements(xc)-1 DO BEGIN
         IF mask[xc(j),yc(j)] eq 1. THEN BEGIN  
            dw_temp=reform(data[xc(j),yc(j),0:40]) 
            IF n_elements(phold) EQ 0 THEN phold=dw_temp ELSE phold=[[phold],[dw_temp]]
         ENDIF 
     ENDFOR
     
     If n_elements(phold) GT 0 THEN BEGIN
         temp.phi=phi_ang*(i-1)
         temp.mn=mean(phold)
         temp.med=median(phold)
         temp.sd=sqrt((moment(phold))[1])
        ; plothist,phold,bin=1,xhist,yhist
        ; plots,[1,1]*temp.mn,[0,max(yhist)],color=0
        ; plots,[1,1]*temp.med,[0,max(yhist)],color=0
         

         IF n_elements(results_dw) EQ 0 THEN results_dw=temp ELSE results_dw=[temporary(results_dw),temp]
    ENDIF
     
     counter,i,no_phi, 'Calculating spectra '

ENDFOR
save,results_dw,filename=outpath+'aa_dw_results.idlsav'



END