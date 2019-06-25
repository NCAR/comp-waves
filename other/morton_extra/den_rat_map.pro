;PURPOSE: - Calculate maps of density and line ratio averaged over a number of files
;
;INPUTS:- date  - date of aligned cube_ivw
;
;OPTIONAL INPUTS: - maxoffset - maximum offset from cross-correlation in wave_tracking.pro
;                 - trange - 2 element vector containing first and last image to use for
;                            calculation
;
;OUTPUTS: - den_map - map of density values
;           lr_map  - map of line ratio values
;
;CALLS: frad_coord.pro, calc_den.pro
;



;----------------------------------------------------------------------------
PRO den_rat_map,date,index=index,den_map=den_map,lr_map=lr_map,$           
                            rplot=rplot,nfiles=nfiles,trange=trange


;##############################################
; LOAD FILES
;##############################################

outpath = 'data/CoMP/wave_tracking_output/'+date+'/'
inpath  = 'data/CoMP/'+strmid(date,8,4,/reverse)+'/'+strmid(date,3,2,/reverse)+'/'+strmid(date,1,2,/reverse)+'/'+date+'.comp.1074.daily_dynamics.5/'
files=find_files('*.fts.gz',inpath)

nfiles=100
test=readfits(files(0),ext=1,/silent)
sz=size(test)
data1074=fltarr(sz(1),sz(2),nfiles)
data1074[*,*,0]=test
FOR i=1,nfiles-1 DO data1074[*,*,i]=readfits(files(i),ext=1,/silent)

inpath  = 'data/CoMP/'+strmid(date,8,4,/reverse)+'/'+strmid(date,3,2,/reverse)+'/'+strmid(date,1,2,/reverse)+'/'+date+'.comp.1079.daily_dynamics.5/'
files=find_files('*.fts.gz',inpath)

nfiles2=40
test=readfits(files(0),ext=1,/silent)
sz=size(test)
data1079=fltarr(sz(1),sz(2),nfiles2)
data1079[*,*,0]=test
FOR i=1,nfiles2-1 DO data1079[*,*,i]=readfits(files(i),ext=1,/silent)

mreadfits,files(0),index

;##############################################
;Set key parameters
;##############################################

nx=sz(1)
ny=sz(2)
xscale=index(0).cdelt1 ;pixel size in arcsec
CASE n_elements(trange) OF
     0: BEGIN
        t1=0 & t2=nfiles-1
        END
     2:BEGIN
        t1=trange[0] & t2=trange[1]
        END
ENDCASE
; Upper and lower limit for the radius in
; !be aware that the maximum offset of the cross-correlation
; (if performed) is added to this value!
lower_r = 228.*xscale
upper_r = 288.*xscale

IF n_elements(maxoffset) LT 1 THEN maxoffset=0. ELSE maxoffset=maxoffset*xscale



image=data1079[*,*,0]

;Calculate maps or r and phi values
r=frad_coord(index(0),image,phi=phi)

; create mask
; mask out pixels to invert
mask=intarr(nx,ny)

good=where(r ge lower_r+maxoffset and r le upper_r)
mask(good)=1.

data_mask1074=data1074
data_mask1079=data1079
for i=0,nfiles-1 do data_mask1074[*,*,i]=data1074[*,*,i]*mask
for i=0,nfiles2-1 do data_mask1079[*,*,i]=data1079[*,*,i]*mask


av_1079=rebin(data_mask1079[*,*,0:19],nx,ny,10)
av_1074=rebin(data_mask1074[*,*,0:49],nx,ny,10)
lr_cube=av_1079/av_1074


lr_map=fltarr(nx,ny,10)
den_map=fltarr(nx,ny,10)

no_r=floor((upper_r-lower_r-maxoffset)/xscale)
rplot=findgen(no_r)*xscale+lower_r+maxoffset


;File that contains values of: 
;line ratio vs height - ratio_array
;height -h 
;density - den
;
restore,'data/comp/line_ratio/rat.sav'

;Calculate density & line ratio for each point

FOR k=0,9 DO BEGIN
FOR i=1,no_r-1 DO BEGIN

   ring=where(r gt (lower_r+maxoffset+(i-1)*xscale) and r lt (lower_r+maxoffset+i*xscale) )
   xc=ring mod nx
   yc=ring/nx

 
  FOR j=0,n_elements(xc)-1 DO BEGIN
  
        ratio=lr_cube[xc(j),yc(j),k] 
        dense=calc_den(ratio_array,h,den,ratio,rplot(i))
        lr_map[xc(j),yc(j),k]=ratio
        den_map[xc(j),yc(j),k]=dense  
   
        
        
     ;ENDIF 
     
  ENDFOR
ENDFOR
ENDFOR

set_plot,'ps'


loadct,4,/silent
device,/encapsul,/color,filename=outpath+'lr_map.eps'
tvim,lr_map[*,*,0]>0.15<0.6
cgcolorbar,range=[0.15,0.7]
device,/close

device,/encapsul,/color,filename=outpath+'den_map.eps'
tvim,den_map[*,*,0]>7.2<8.5
cgcolorbar,range=[7.2,9.2]

device,/close

loadct,0,/silent
!p.multi=0
set_plot,'x'


END