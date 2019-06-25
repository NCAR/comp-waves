;+
; :Description:

;   procedure to compute analytic fit to a gaussian sampled at three points

;INPUTS:  i1, i2, i3 -  the intensity images/cubes increasing monotonically in wavelength
;         d_lambda -  the wavelegth spacing of the samples
;
;OUTPUTS: doppler_shift - the shift of the fit gaussian from the center wavelength in the 
;                         same units as d_lambda
;         width -  the linewidth in the same units as d_lambda
;         i_cent - the central intensity of the gaussian in the same units as i1, i2, i3
;
;OPTIONAL INPUTS: /get_rid - remove NAN's
;
; :Params:
;    profile
;    d_lambda
;    doppler_shift
;    width
;    i_cent
;
;
;
;HISTORY:  Analytic_gauss_fit written by Christian Bethge
;          07/2014 RJM updated to use arrays rather than individual elements -x200 faster
;                       
;
;-
pro analytic_gauss_fit_fast,I1,I2,I3,d_lambda,doppler_shift,width,i_cent,get_rid=get_rid

  sz=size(I1)
  nx=sz(1)
  ny=sz(2)

  ;FOR 2d or 3d arrays
  IF sz(0) EQ 3 THEN BEGIN
     nt=sz(3)
     doppler_shift=dblarr(nx,ny,nt) & width=dblarr(nx,ny,nt) & i_cent=dblarr(nx,ny,nt)
  ENDIF ELSE BEGIN
     doppler_shift=dblarr(nx,ny) & width=dblarr(nx,ny) & i_cent=dblarr(nx,ny)
  ENDELSE


  a=alog(I3/I2) & b=alog(I1/I2)


  width=sqrt( -2D*d_lambda^2D/(a+b) )
  doppler_shift=width^2D/(4D*d_lambda)*(a-b)
  i_cent=i2*exp(doppler_shift^2D/width^2D)

  
  width[where(I1 lt 0. or I2 lt 0. or I3 lt 0.)]=0D
  doppler_shift[where(I1 lt 0. or I2 lt 0. or I3 lt 0.)]=0D
  i_cent[where(I1 lt 0. or I2 lt 0. or I3 lt 0.)]=0D 

  
  
  width[where(-2D*d_lambda^2D/(a+b) lt 0)]= 0D 
  doppler_shift[where(-2D*d_lambda^2D/(a+b) lt 0)]= 0D
  i_cent[where(-2D*d_lambda^2D/(a+b) lt 0)]= 0D 

  IF keyword_set(get_rid) THEN BEGIN
     width[where(finite(width,/nan,sign=1))]=0D 
     width[where(finite(width,/nan,sign=-1))]=0D
     
     doppler_shift[where(finite(doppler_shift,/nan,sign=1))]=0D 
     doppler_shift[where(finite(doppler_shift,/nan,sign=-1))]=0D

     i_cent[where(finite(i_cent,/nan,sign=1))]=0D 
     i_cent[where(finite(i_cent,/nan,sign=-1))]=0D
     i_cent[where(finite(i_cent,/infinity))]=0D
  
  ENDIF
  
  
end
