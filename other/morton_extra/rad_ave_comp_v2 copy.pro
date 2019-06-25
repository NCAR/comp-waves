;PURPOSE: - Calculate radially and angular-ly averaged power spectra for CoMP data
;
;INPUTS:- index - initial of observations (only used to calculate the radial coordinates) 
;         date  - date of aligned cube_ivw
;
;OPTIONAL INPUTS: - maxoffset - maximum offset from cross-correlation in wave_tracking.pro
;                 - apod - degree of apodisation (<1.) - default is 0.1
;
;                   extra - string relating to additional characters in cube_ivw name
;                   /dopwid - calculates average using Doppler width
;                   /inten - calculates average using intensity 
;                   /errors - will find error file for the day and plot power spectra for the errors
;                   /do_phi - will provide the data averaged over different phi angles (no radial seperation)
;
;
;OUTPUTS: - ravpow - radially averaged power array of frequency vs height
;           phipow - angular-ly averaged power array of frequency vs angle
;           rplot  - array of radial coordinates in arcsec
;
;TO DO: - 
;         Problem with 'duff' time-series - current removal of time-series with number of zeros is problematic
;         Might be better to exclude series on rms power values
;



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


 ; apodizes time-series, with optional de-trending 
;----------------------------------------------------------------------------
FUNCTION apod3dcube,cube,apod,mask,cpg=cpg ;,detrend=detrend
    
  ; get cube dimensions
  sizecube=size(cube)
  nx=sizecube[1]
  ny=sizecube[2]
  nt=sizecube[3]
  apocube = fltarr(nx,ny,nt)

  ; define temporal apodization
  apodt = fltarr(nt)+1
  IF (apod ne 0) THEN BEGIN
    apodrimt = nt*apod
    apodt[0] = (sin(!pi/2.*findgen(apodrimt)/apodrimt))^2
    apodt = apodt*shift(rotate(apodt,2),1)   
  ENDIF
  CPG=total(apodt)/n_elements(apodt) ;Coherent Power Gain - correction factor needed for power after apodisation

  ;temporal mean removal and apodization
  apocube=(cube-rebin(mean(cube,dim=3),nx,ny,nt))*transpose(rebin(apodt,nt,nx,ny),[1,2,0])

  return,apocube
END

;Taken fro average_spec.pro and changed
;may cause problems if both are run together
;

FUNCTION power_fit_cube,nt,input,noise=noise


tspec=fltarr(nt/2,3)
nelm=n_elements(input[0,*])

IF n_elements(noise) EQ 1 THEN start=[200,-10.,1,0] ELSE start=[200,-6.,1,0]

;Forget about DC component
FOR i=1,nt/2-1 DO BEGIN
    
    plothist,alog(input[i,*]),xhist,yhist,/noplot
    IF (where(yhist eq 0))[0] NE -1 THEN yhist(where(yhist eq 0))=1    
    reg=mpfitfun('mygauss',xhist,yhist,sqrt(yhist),start,/quiet,perror=pp,bestnorm=bn,dof=dof)
    a=moment(alog(input[i,*])) ; calculates moments for log data - shows    log-normal behaviour
    
    tspec[i,0]=reg(1); a(0) ;exp(a[0]) 
    tspec[i,1]=(pp[1]);*sqrt(bn/dof) ;a[1]
    tspec[i,2]=reg(2);sqrt( (exp(tspec[i,1]^2)-1)*exp(2*tspec[i,0]+tspec[i,1]^2))

ENDFOR
return,tspec

END


;----------------------------------------------------------------------------
PRO rad_ave_comp_v2,date,index=index,ravpow=ravpow,phiavpow=phiavpow,$           
                            rplot=rplot,apod=apod,extra=extra,dopwid=dopwid,inten=inten,errors=errors,$
                            do_phi=do_phi,do_atrous=do_atrous,phi_ang=phi_ang


; paths where the CoMP data is located and
; where the .sav-files should be written to
inpath  = 'analysis/CoMP/wave_tracking_output/'+date+'/'
outpath = 'analysis/CoMP/wave_tracking_output/'+date+'/'

;inpath_ind  = 'data/CoMP/'+strmid(date,8,4,/reverse)+'/'+strmid(date,3,2,/reverse)+'/'+strmid(date,1,2,/reverse)+'/'
;files=find_files('*.fts.gz',inpath_ind)
;mreadfits,files(0),index


IF keyword_set(errors) THEN restore,inpath+'cube_ivw_errors.sav',/verbose ELSE $
IF NOT keyword_set(extra) THEN restore,inpath+'cube_ivw_'+date+'.sav',/verbose $
ELSE restore,inpath+'cube_ivw_'+date+'_'+extra+'.sav',/verbose

IF keyword_set(errors) THEN BEGIN
    data=dv<100.
    IF keyword_set(dopwid) THEN data=dw
    IF keyword_set(inten) THEN data=di
    

ENDIF ELSE BEGIN 
    data=cube_v
    IF keyword_set(dopwid) THEN data=cube_w
    IF keyword_set(inten) THEN data=cube_i
ENDELSE


sz=size(data)
nx=sz(1)
ny=sz(2)
nt=sz(3)
xscale=index.cdelt1 ;pixel size in arcsec


; Upper and lower limit for the radius in
; !be aware that the maximum offset of the cross-correlation
; (if performed) is added to this value!
lower_r = index.lower_r*xscale
upper_r = index.upper_r*xscale
maxoffset=index.maxoffset*xscale


r=frad_coord(index,data[*,*,0],phi=phi)

; create mask
; mask out pixels to invert
mask=index.mask
data_mask=data
for i=0,nt-1 do data_mask[*,*,i]=data[*,*,i]*mask

no_r=floor((upper_r-lower_r-maxoffset)/xscale)
rplot=findgen(no_r)*xscale+lower_r+maxoffset

ravpow=fltarr(nt/2,no_r,2)
rpow_elem=0.

print,'Starting'

;Apodise the time-series for each pixel
IF n_elements(apod) LT 1 THEN apod=0.1
apocube=apod3dcube(data_mask,apod,mask,cpg=cpg)


;Very high pass filter
;Good at isolating noise and keeping small scales features
;
IF keyword_set(do_atrous) THEN BEGIN
   FOR i=0,nt-1 DO BEGIN
       atrous,apocube[*,*,i],decomposition=a,n_scales=1
       oas=a[*,*,1]-smooth(a[*,*,1],3)
       ;oas=oas-smooth(oas,3)
       apocube[*,*,i]=temporary(apocube[*,*,i])-oas
   ENDFOR
ENDIF



IF NOT keyword_set(do_phi) THEN BEGIN
        ;Radial averaging
        FOR i=1,no_r DO BEGIN

           ring=where(r gt (lower_r+maxoffset+(i-1)*xscale) and r lt (lower_r+maxoffset+i*xscale) )
         ;print,lower_r+maxoffset+(i-1)*xscale,lower_r+maxoffset+(i)*xscale
         ;print,n_elements(ring)
          xc=ring mod nx
          yc=ring/nx

          phold=fltarr(nt/2)
          FOR j=0,n_elements(xc)-1 DO BEGIN
          
             no_zero=n_elements(where(apocube[xc(j),yc(j),*] eq 0.))
              
             IF no_zero lt 5 THEN BEGIN 
                times=apocube[xc(j),yc(j),*]
                

                ;Excludes pixels with time-series that have spikes greater
                ;than pm 4 sigma from the mean
                
                vals=moment(times)
                cl=(vals[0]+4.*sqrt(vals[1]))
                cll=(vals[0]-4.*sqrt(vals[1]))
                IF n_elements(where(times GT cl OR times LT cll)) GT 1 THEN BEGIN
                  
                ENDIF ELSE BEGIN 
                   
                       dft=fft(times,-1) 
                       dft=temporary(dft[0:nt/2-1])
                       phold=[[phold],[2.*abs(dft)^2]]

                ENDELSE

             ENDIF 
             
          ENDFOR
          rpow_elem=[[rpow_elem],[n_elements(phold[0,*])]]
          pelm=n_elements(phold[0,*])
           
         
          FOR jj=0,nt/2-1 DO BEGIN
              med=median(phold[jj,1:pelm-1])
              sd=sqrt((moment(phold[jj,1:pelm-1]))[1])
              ravpow[jj,i-1,0]=med
              ravpow[jj,i-1,1]=sd
              ;print,jj,i-1
          ENDFOR
         
        ENDFOR




ENDIF ELSE BEGIN

         ;Phi averaging
          IF n_elements(phi_ang) EQ 0 THEN phi_ang=45
          no_phi=360/phi_ang
          phiavpow=fltarr(nt/2,no_phi)
          phipow_elem=0. 

          FOR i=1,no_phi DO BEGIN

               
               wedge=where(phi gt phi_ang*(i-1) and phi lt phi_ang*i and mask gt 0)
               xc=wedge mod nx
               yc=wedge/nx
               
               phold=fltarr(nt/2)
               FOR j=0,n_elements(xc)-1 DO BEGIN
            
                   no_zero=n_elements(where(apocube[xc(j),yc(j),*] eq 0.))
                    
                   IF no_zero lt 5 THEN BEGIN 
                      dft=fft(apocube[xc(j),yc(j),*],-1)  
                      dft=temporary(dft[0:nt/2-1]) 
                     
                      phold=[[phold],[2.*abs(dft)^2/cpg^2]]
                     
                   ENDIF 
               
               ENDFOR

            res=power_fit_cube(nt,phold,noise=noise)  
            ;bootstrap means
            num_mc=500
            num_bs=n_elements(phold[0,*])
            bs_arr=round(randomu(systime,num_bs,num_mc)*(num_bs-1))
            FOR kk=0,num_mc-1 DO BEGIN
                IF n_elements(mn_bs) EQ 0 THEN mn_bs=mean(phold[*,bs_arr[*,kk]],dim=2)$
                ELSE mn_bs=[[mn_bs],[mean(phold[*,bs_arr[*,kk]],dim=2)]]

            ENDFOR
   ;         phiavpow[*,i-1]=total(phold[*,1:n_elements(phold[0,*])-1],2)/n_elements(phold[0,*])
            stop

          ENDFOR
   

          ;Plotting of image and angle contours

          ;Define colour table
          aia_lct,wavelnth='193',/load,/silent
          tvlct,r,g,b,/get
          cols=[[r],[g],[b]]

          yran=([0,ny-1]-index.crpix2)*xscale
          xran=([0,nx-1]-index.crpix1)*xscale
          intplot=image(reform(cube_i[*,*,0])*mask,axis_style=0,rgb_table=cols)
          xaxis=axis('x',location=0,coord_trans=[xran[0],xscale],title='Distance (arcsec)')
          xaxis=axis('y',location=0,coord_trans=[yran[0],xscale],title='Distance (arcsec)')

          cvalue=phi_ang*findgen(no_phi)
          contplot=contour(phi,overplot=intplot,c_value=cvalue,c_linestyle=2,color='white',c_label_show=0)
          
          xloc=index.crpix1+(280.)*sin(cvalue*!const.dtor)
          yloc=index.crpix2+280.*cos(cvalue*!const.dtor) 
          mylabels=text(xloc,yloc,strtrim(uint(cvalue),2),overplot=contplot,color='white',/data)
          contplot.Scale,0.9,0.9
          contplot.Save,'test.eps'

ENDELSE
stop

END