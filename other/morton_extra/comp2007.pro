
;goto,ratiomap
;;;;; PRO is in Ret is out 
restore,'analysis/comp/2007/cube_20051030.sav',/v
restore,'analysis/comp/2007/waveangle_2005.idlsav',/v

sz=size(cube_v,/dim)
nx=sz(0)
ny=sz(1)
nt=sz(2)
ns=4.

x=rebin(findgen(nx)-nx/2,nx,ny)
y=transpose(x)
rad=sqrt(x^2+y^2)




dx      = 20      ;
nbox    = 2*dx+1  ;box size in pixels
nsmooth = 15 ;?15     ;smoothing width for cross spectra (nominally 11)
limit   = 0.5     ;set coherence limit for acceptance
npix    = 10      ;minimum number of coherent pixels
fwidth = 0.0015
freq_filt=0.0035

nspec = nt/2
spec=complexarr(nx,ny,nspec)
;tp=time
cad=29.
wave_angle=fltarr(nx,ny)
angle_error=fltarr(nx,ny)
angle_error2=fltarr(nx,ny) ; For testing RJM
coh_measure = fltarr(nx,ny,2)
f=8

nspec=nt/2

freq   = findgen(nspec)/(float(nspec*2)*cad)
filter = exp( -(freq-freq_filt)^2/fwidth^2 )
filter(0) = 0.                      ; set dc to zero
filter(where (freq lt .001)) = 0.   ; set low frequencies to zero
filter=filter/total(filter)


;reduced for 2005 cube
ny1=400
nx1=300




;restore,'analysis/CoMP/wave_tracking_output/20101027/cube_ivw_20101027.sav'
;restore,'analysis/CoMP/wave_tracking_output/20101027/wave_angle_20101027_3_3.5mHz.sav'
loadct,13

;goto, ratiomap
;read,try,prompt='enter try'

window,3,xs=1000,ys=700
;velocity=cube_v
plot_image,cube_v(*,*,0)



;stop
;; calculating the radius from 3 points;;;;;;

;print,'pick points on occulting disk edge'

pick,x1,y1
pick,x2,y2
pick,x3,y3
xn=[x1,x2,x3]
yn=[y1,y2,y3]
cir_3pnt,xn,yn,r0,x0,y0
;restore,'center.sav',/v
print,'radius is ',r0 
;print,'center is', x0,y0


;;;;;;;;;;;;;;;;;;;;;;

;;;;;; defining mask for cube_vimages;;;;;


lower_r=r0+5
upper_r=lower_r+80.

good=where(rad ge lower_r and rad le upper_r)

;
;
;mask1=index.mask
mask1=fltarr(nx,ny)
mask1[good]=1.

S = REPLICATE(1, 3, 3) ;;; image operator


;mask1=dilate(mask1,S)


plot_image,mask1
print,'mask done!'

goto, partspeed

loadct, 13, /silent


print,'Begin FFT of data cube'
;goto,partspeed




h=do_apod(nt,cpg)

FOR iy=0,ny1-1 DO FOR ix=0,nx1-1 DO BEGIN IF mask1(ix,iy) EQ 1 THEN BEGIN

  d=reform(cube_v[ix,iy,*]-mean(cube_v[ix,iy,*])) 
  sp=fft(d*h)    
  spec[ix,iy,*]=sp[0:nspec-1]
endif
;save,spec,filename='spec1.sav'
ENDFOR


xbox=rebin(findgen(nbox)-dx,nbox,nbox)  ;coordinates of points in box
ybox=transpose(xbox)

pixcounter = 0L
old_perc_proc = 0
tot_pix = long(where(mask eq 1))
conjspec=conj(spec)
lar_g=real_part(smooth(spec*conjspec,[1,1,nsmooth]))

print,'begin angle calc'

FOR iy=0,ny1-1 DO $ 
FOR ix=0,nx1-1 DO BEGIN
  IF mask1[ix,iy] EQ 1 THEN BEGIN
  ;================================================================
  ;  compute coherence for contiguous pixels with coherence greater than limit
  ;================================================================

  coh=fltarr(nbox,nbox)           ;initialize coherence array to zero
  coh_mask=intarr(nbox,nbox)    ;initialize mask of good points in coherence
  coh(nbox/2,nbox/2)=1.
  coh_mask(nbox/2,nbox/2)=1.  ;;remove coherence mask 

  spec1=reform(spec[ix,iy,*])
  g1=reform(lar_g[ix,iy,*])

  count=1
  ncoh=1
  WHILE count GT 0 AND ncoh LT 100 DO BEGIN
        cont=where(coh_mask eq 1.)
        cnew=[cont-1,cont+1,cont+nbox,cont-nbox]  ;test pixels adjacent to contiguous pixels

        new_mask=intarr(nbox,nbox)+1   ;identify new contiguous pixels
        new_mask(cnew)=1.
        new_mask=new_mask*(1.-coh_mask)

        new=where(new_mask gt 0.,new_count)
       
        count=0
        FOR ic=0,new_count-1 DO BEGIN
            icx=new(ic) mod nbox
            icy=new(ic)/nbox

            ii=ix-nbox/2+icx
            jj=iy-nbox/2+icy

            IF ii LT 0 OR ii GT nx-1 OR jj LT 0 OR jj GT ny-1 THEN coh(icx,icy)=0. ELSE BEGIN
        	  IF coh(icx,icy) EQ 0. THEN BEGIN
          	         spec2c=reform(conjspec(ii,jj,*))
                     cspec=spec1*spec2c
                     cspec=smooth(cspec,nsmooth)
                     g2=reform(lar_g[ii,jj,*])
                                                               
                     coh(icx,icy)=mask1(ii,jj)*total( filter*( abs(cspec)/sqrt(g1*g2) ) ) ;
                     ENDIF

                  IF coh(icx,icy) GT limit THEN BEGIN
          	     coh_mask(icx,icy)=1.
          	     count=count+1
                  ENDIF
            ENDELSE
         ENDFOR

         good=where(coh_mask eq 1.,ncoh)

   ENDWHILE

   IF (ncoh gt npix) THEN BEGIN
     island = coh*coh_mask
     ellipsepts = fit_ellipse(where(island ne 0.),axes=axes,xsize=nbox,ysize=nbox)
     coh_measure[ix,iy,*] = axes
     
   ENDIF

;=============================================================== 
;  find angle of coherence island if enough good points  
;===============================================================

   IF ncoh GT npix THEN BEGIN

 
       weight=coh(good)^2
       theta=min_dist_fit(xbox(good),ybox(good),error=error) ; analytical solution to the minimum distance problem

       IF keyword_set(altangle) THEN BEGIN

           xe=xbox(good) & ye=ybox(good)
           sxy=total((xe-mean(xe))*(ye-mean(ye)))

           IF round(sxy) NE 0. THEN BEGIN
    
                 sixlin,xbox(good),ybox(good),a,siga,b,sigb
                 theta=[atan(b[2])/(!dtor)]
                 error=[atan(sigb[2])/(!dtor)]
           ENDIF
       ENDIF
        
       wave_angle(ix,iy)=theta
       angle_error(ix,iy)=error

    ENDIF
ENDIF
ENDFOR

;save,wave_angle,mask1,filename='/analysis/comp/2007/waveangle_2007'+string(npt_init,format='(I4.4)')+'.sav'
;set_plot,'PS'

;device,filename='steve_plot1.eps',/color,bits_per_pixel=16,xsize=75,ysize=25

;plot,xbox(good),ybox(good),psym=6,xr=[-dx,dx],yr=[-dx,dx]
;oplot,xfit,yfit
;device,/close

stop

;set_plot,'PS'

;device,filename='steve_plot2.eps',/color,bits_per_pixel=16,xsize=75,ysize=25
;tv,bytscl(rebin(coh*coh_mask,f*nbox,f*nbox,/sample),0.,1.)


;device,/close

set_plot,'X'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

partspeed:

angle=wave_angle
xscale= 3.2267284

norm_cadence =cad

;###############################################
;This makes cube odd length - easier to split FFT
;as it gets rid of Nyquist vapue
if ((nt mod 2) eq 0) then nt = nt-1 else nt = nt
cube_v=temporary(cube_v[0:nx-1,0:ny-1,0:nt-1])

;;; *** choose length of track here ***
npt_init=61 ;number of steps to map along track (make odd) 
dx=1.0      ;step size along track

;===========================
npt= npt_init
ang_track=fltarr(npt)   ;angle along track
vmap=fltarr(nt,npt)     ;array to hold space-time diagram

f=findgen(nt)/(float(nt)*norm_cadence)       ;compute temporal frequency scale
for i=0,nt/2 do f(nt-i-1)=f(i+1)
freq=rebin(f,nt,npt)                         ;along track,

sf=findgen(npt)/(float(npt)*xscale)          ;compute spatial frequency scale ;length is in Mm (res=4.46*0.725 Mm)
for i=0,npt/2 do sf(npt-i-1)=sf(i+1)         ;make retrograde frequencies negative
sp_freq=rebin(transpose(sf),nt,npt)          ;along track

phase_speed=freq/sp_freq                     ;compute phase speed = (w/k)



w=0.001
filter=1.-exp(-freq^2/w^2)

angle=median(angle,3)   ;;removed for various tracks
pro_pow_f=fltarr(nx,ny,nt/2)
ret_pow_f=fltarr(nx,ny,nt/2)
pow_rat=fltarr(nx,ny,nt/2)

sign_power=fltarr(nx,ny)
pro_speed=fltarr(nx,ny)
ret_speed=fltarr(nx,ny)
phase_speed=fltarr(nx,ny)
pro_speed_err=fltarr(nx,ny)
ret_speed_err=fltarr(nx,ny)
phase_speed_err=fltarr(nx,ny)
pro_no_tr=fltarr(nx,ny)
ret_no_tr=fltarr(nx,ny)
ystart=1
xstart=1
yend=ny1-1
xend=nx1-1
tot_pix=(xend-xstart)*1L*(yend-ystart)

;TIC

FOR iy=ystart,yend DO BEGIN
  FOR ix=xstart,xend DO BEGIN
   ; clock = TIC('FFT') & $
     IF ix le nx/2-1 THEN hemisphere=1 ELSE hemisphere=-1
     
     IF mask1[ix,iy] EQ 1 THEN BEGIN
         
            img1 = cube_v[0:nx-1,0:ny-1,0]   ;reinitialize velocity image
            img1 = bytscl(img1,-8,8,top=254)
            img2 = bytscl(angle,-90,90,top=254)
             
             npt=npt_init         ;every times we have to do this
             xtrack=fltarr(npt)   ;x and y position along track
             ytrack=fltarr(npt)   ;............................
             ang_track=fltarr(npt)   ;angle along track

             imid=fix(npt/2)

             ;##############################################################
           ; use angle at each point to define track, starting at cursor position as central point
           ;
           ;Place initial position +angle at centre of arrays xtrack + ytrack +ang_track
           ;##############################################################
    
           xtrack(imid)=float(ix) & ytrack(imid)=float(iy)
           ang_track(imid)=angle(ix,iy)

           ntrack=0        ;this index will count the number of pixels on the track

         ;this condition gives the lower limit of the track
           IF sqrt( (float(ix)-nx/2.+0.5)^2 + (float(iy)-ny/2.+0.5)^2 )-imid LT lower_r THEN $
                   short_track=1 else short_track=-1 
  
         ;##############################################################
         ;  first, move out from cursor position
         ;##############################################################
      
         FOR i=imid+1,npt-1 DO BEGIN

               ;Calculates which x + y direction to follow
               ;Uses angle calculated from cross-correlation in wave_tracking
               ;Use of hemisphere ensure signal always steps away from occulting disk
               ;This means value should always greater than 0 and less than nx-1 (ny-1)
               xtrack(i)=0. > xtrack(i-1)-hemisphere*dx*cos( ang_track(i-1)*!dtor) < float(nx-1)
               ytrack(i)=0. > ytrack(i-1)-hemisphere*dx*sin( ang_track(i-1)*!dtor) < float(ny-1)

               intx=0 > fix(xtrack(i)+0.5) < (nx-1)
               inty=0 > fix(ytrack(i)+0.5) < (ny-1)
            

                     ;this means if you are near the lower boundary
                     pixcheck=sqrt((intx-nx/2.+0.5)^2 + (inty-ny/2.+0.5)^2)
               IF pixcheck lt lower_r THEN BEGIN
                        ;this will not allow to track further down in the inward direction!
                  xtrack(i)=0.
                  ytrack(i)=0.

                  BREAK
                        
               ENDIF
           
               IF pixcheck GT upper_r THEN BEGIN
                  xtrack(i)=0.
                  ytrack(i)=0.

                  BREAK
               ENDIF

                 ;IF mask[intx,inty] LT 1 THEN BEGIN ;don't use bad pixels
               ;   xtrack(i)=0.
               ;   ytrack(i)=0.
               ;   BREAK
               ;ENDIF

               ntrack=ntrack+1

                 IF keyword_set(debug) THEN BEGIN
              img1(intx,inty)=255 ;put track into displayed images
              img2(intx,inty)=255
               ENDIF
               
                 ang_track(i)=angle(intx,inty)  ; nearest neighbor interp for angle
                                           ; (works better than bilinear at boundaries)

                 ; if angle is zero (not defined), use previous value
               ;IF ang_track(i) eq 0. THEN ang_track(i)=ang_track(i-1)

               ; if big difference in adjacent angles, resolve 180 degree ambiguity
               ; by choosing angle closest to ang_track(i-1)
               IF abs(ang_track(i)-ang_track(i-1)) gt 90. THEN BEGIN
                 IF ang_track(i)-ang_track(i-1) GT 0. THEN ang_track(i)=ang_track(i)-180. ELSE ang_track(i)=ang_track(i)+180.
               ENDIF

         ENDFOR

         ;##############################################################
         ;  next, move in from cursor position
         ;##############################################################
         FOR i=imid-1,0,-1 DO BEGIN
              factor=(ang_track(i+1)+180.)*!dtor
              xtrack(i)=0. > xtrack(i+1)-hemisphere*dx*cos(factor ) < float(nx-1)
              ytrack(i)=0. > ytrack(i+1)-hemisphere*dx*sin( factor) < float(ny-1)

              intx=0 > fix(xtrack(i)+0.5) < (nx-1)
              inty=0 > fix(ytrack(i)+0.5) < (ny-1)

              ;this means if you are near the lower boundary
              pixcheck=sqrt((intx-nx/2.+0.5)^2 + (inty-ny/2.+0.5)^2)
              IF pixcheck lt lower_r THEN BEGIN
                 ;this will not allow to track further down in the inward direction!
                xtrack(i)=0.
                ytrack(i)=0.
              BREAK ; stops tracks when hits lower boundary
                ENDIF
           
              IF pixcheck GT upper_r THEN BEGIN
                  xtrack(i)=0.
                  ytrack(i)=0.
                BREAK ; stops tracks when hits upper boundary
              ENDIF

                   ;  IF mask[intx,inty] LT 1 THEN BEGIN ;don't use bad pixels
             ;     xtrack(i)=0.
             ;     ytrack(i)=0.
             ;     BREAK
             ;  ENDIF

              ntrack=ntrack+1

              ang_track(i)=angle(intx,inty)  ; nearest neighbor interp for angle
                                           ; (works better than bilinear at boundaries)

              ; if angle is zero (not defined), use previous value
              ;IF ang_track(i) EQ 0. THEN ang_track(i)=ang_track(i+1)

              ; if big difference in adjacent angles, resolve 180 degree ambiguity
                    ; by choosing angle closest to ang_track(i-1)
              IF abs(ang_track(i)-ang_track(i+1)) GT 90. THEN BEGIN
                IF ang_track(i)-ang_track(i+1) GT 0. THEN ang_track(i)=ang_track(i)-180. ELSE ang_track(i)=ang_track(i)+180.
              ENDIF

         ENDFOR


         IF ntrack LT 4 THEN CONTINUE ;this tells that we are taking at least 5 points in the track
         
         IF ntrack NE ntrack/2*2+1 THEN ntrack=ntrack/2*2+1  ;Gives the total number of pixels in the track.
                                                             ;For normal case it is npt, but when the track
                                                             ;is shorter, it is < npt
                                 ;This also makes npt odd
         IF ntrack LT npt THEN BEGIN    ;......If the track is shorter, then it will do the followings
           in=where(xtrack ne 0)
           newxtrack=xtrack(in) ;......it creates new track of shorter length
           xtrack=newxtrack
           newytrack=ytrack(in)
           ytrack=newytrack
           npt=n_elements(in)
               ang_track=ang_track(in)
         ENDIF ELSE BEGIN
                 npt=npt_init                    ; We need the following lines for the new track
         ENDELSE


         vmap=fltarr(nt,npt)           ;Again we need to make it (array to hold space-time diagram)

         ; interpolate velocity onto track
         FOR i=0,nt-1 DO vmap[i,0:npt-1]=interpolate(cube_v[0:nx-1,0:ny-1,i],xtrack,ytrack,cubic=-0.5)

         vmap=vmap-mean(vmap)
         vmap=vmap-transpose(rebin(mean(vmap,dim=1),npt,nt),[1,0]) ;remove temporal mean
   

         trans=fft(vmap,-1)       ;compute fourier transform
         

         pro_trans=trans             ;select prograde waves (assume nt even, npt odd)
         pro_trans(1:nt/2-1,0:npt/2)=0.
         pro_trans(nt/2+1:nt-1,npt/2+1:npt-1)=0.
         pro_vel=real_part(fft(pro_trans,1))
    	   sd=where(ytrack eq iy and xtrack eq ix)
    	   provelt=pro_vel(*,sd) ;;;; check frequency + /-
    	   provelfft=fft(provelt,-1,/center) ;;; -1 vs 1 fft idl doc
    	
         pro_pow_f(ix,iy,*)=abs(provelfft[0:nt/2-1])^2
         

         ret_trans=trans             ;select retrograde waves
         ret_trans(1:nt/2-1,npt/2+1:npt-1)=0.
         ret_trans(nt/2+1:nt-1,0:npt/2)=0.
         ret_vel=real_part(fft(ret_trans,1))
    	   ;sd1=size(ret_vel,/dim)
    	   retvelt=ret_vel(*,sd) ;;;; check freq values 
    	   retvelfft=fft(retvelt,-1,/center)
    	
         ret_pow_f(ix,iy,*)=abs(retvelfft[0:nt/2-1])^2
         
         
          
          ;r=compute_speed(pro_vel,xscale,norm_cadence)
          ;pro_speed(ix,iy)=r(0)
          ;pro_speed_err(ix,iy)=r(1)
          ;pro_no_tr[ix,iy]=no_tr
      
          ;r=compute_speed(ret_vel,xscale,norm_cadence,/ret)
          ;r=compute_speed_new(ret_vel,xscale,norm_cadence,/ret,no_tr=no_tr)
          ;ret_speed(ix,iy)=r(0)
          ;ret_speed_err(ix,iy)=r(1)
            ;ret_no_tr[ix,iy]=no_tr

          ;phase_speed(ix,iy)=(pro_speed(ix,iy)/pro_speed_err(ix,iy)^2 + abs(ret_speed(ix,iy))/$
                            ;ret_speed_err(ix,iy)^2)/(1.0/pro_speed_err(ix,iy)^2 + 1.0/ret_speed_err(ix,iy)^2)
          ;phase_speed_err(ix,iy)=(1.0/pro_speed_err(ix,iy)^2 + 1.0/ret_speed_err(ix,iy)^2)^(-0.5)

         ;##############################################################
         ; compute prograde and retrograde power
         ;##############################################################
          ;pro_pow=abs(pro_pw_trans)^2
          ;pro_power(ix,iy)=mean(pro_pow(good_phase))
          ;pro_power(ix,iy)=total(pro_pow)
          ;ret_pow=abs(ret_pw_trans)^2
          ;ret_power(ix,iy)=mean(ret_pow(good_phase))
          ;ret_power(ix,iy)=total(ret_pow)

         ;##############################################################
         ; determine whether track is outward or inward, which determines sign of prograde and retrograde
         ;##############################################################

          sign_power(ix,iy)=-(rad(xtrack(0),ytrack(0))-rad(xtrack(npt-1),ytrack(npt-1)))/ $
                            abs(rad(xtrack(0),ytrack(0))-rad(xtrack(npt-1),ytrack(npt-1)))

          if keyword_set(debug) then !p.multi=0
            
        
     ENDIF
  ENDFOR

   counter,ix-xstart+(xend-xstart)*1L*(iy-ystart),tot_pix,/percent
 ;  toc, clock & $
ENDFOR

;toc

stop

;ratiomap:




pow_rat=fltarr(620,620,512)
temp=fltarr(620,620)
indx=where(mask1 eq 1)
for i=0,512 do begin
    temp1=pro_pow_f[*,*,i]/ret_pow_f[*,*,i]
    temp[indx]=temp1[indx]
    pow_rat[*,*,i]=temp

end
      

;pow_rat(ix,iy,i)=(ret_pow_f(ix,iy,i)-pro_pow_f(ix,iy,i))/(ret_pow_f(ix,iy,i)+pro_pow_f(ix,iy,i))


;ENDFOR
;ENDFOR

;ENDFOR

save,pro_pow_f,ret_pow_f,pro_power,ret_power,sign_power,pro_speed,ret_speed, $
pro_speed_err,ret_speed_err,pow_rat,phase_speed,phase_speed_err,pro_no_tr,ret_no_tr,filename='power_maps2007'+string(npt_init,format='(I4.4)')+'.sav'



;pow_rat=pow_rat(*,*,1:511)

stop





pw=image(-90>wave_angle<90,title='wave_angle',$
          aspect_ratio=0.9,rgb_table=13,margin=0.1,/buffer)
 ax1=axis('x',title='Pixels',location='bottom')  ;,coord_transform=[0:250]*4.46)                                           
 
 ax2=axis('y',title='Pixels',location='left')  ;,coord_transform=[0:350]*4.46)                                             
 
 cb=colorbar(target=pl,rgb_table=13,range=[-90,90],position=[0.5,0.7,0.8,0.72],/normal)
 
pw.save,'wave_angle.png'

END


