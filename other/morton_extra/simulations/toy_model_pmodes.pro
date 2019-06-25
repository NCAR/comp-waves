PRO toy_model_pmodes

outpath='analysis/CoMP/papers/dopp_vel_2/'

files=find_files('*.fts','analysis/mdi_144day/')
mdi_pow=readfits(files)
mdi_freq=findgen(1296)*6.43e-6
mdi_pow=smooth(total(mdi_pow,1),30)

nsamp=500

m=30.
v=55.

phi=sqrt(v + m^2)
mu=alog(m^2/phi)
sigma = sqrt(alog(phi^2/m^2))
x=randomn(1,nsamp)*sigma+mu
y=exp(x)

window,0
plothist,y

gam=5./3.
g=274. ;m/s^2
temp_k=5000. ;Kelvin
h0=!const.k*temp_k/!const.mp/g ;m
cs=sqrt(gam*g*h0) ;m/s
vac=gam*g/4./!pi/cs ;cut-off freq

print,vac

;Radiation modified - Centeno
nfreq=n_elements(mdi_freq)
omega=mdi_freq*2.*!pi
is=fltarr(nfreq)
is[*]=1.
imag=complex(fltarr(nfreq),is)
tau_rad=30. ;s

HR=omega^2*(1+omega^2*tau_rad^2*gam)/g/H0/(1+omega^2*tau_rad^2*gam^2)-1./4/H0^2
HI=(gam-1)*tau_rad*omega^3/(1+omega^2*tau_rad^2*gam^2)/g/H0

kr_sq=0.5*(hr+(hr^2+hi^2)^0.5)
ki_sq=0.5*(-hr+(hr^2+hi^2)^0.5)

;gam_hat=(1.- gam*omega*tau_rad*imag)/(1.-omega*tau_rad*imag)
;c_hat_sq=(gam_hat*g*H0)
;wac_hat_sq=c_hat_sq/4./H0^2
;kz_sq=(omega^2-wac_hat_sq)/c_hat_sq

z=2000e3 ;m - height of chromosphere
amp=mdi_pow*exp(z/2/H0)*exp(-1.*sqrt(ki_sq)*z)


out=fltarr(1296,nsamp,2)
FOR i=0,nsamp-1 DO BEGIN
	ang=cos(y[i]*!dtor)
	mdi_temp2=mdi_pow
    vac_mod=vac*ang
    mdi_temp2[where(mdi_freq lt vac_mod)]=0
    out[*,i,1]=mdi_temp2
    
    h0_mod=!const.k*temp_k/!const.mp/g/ang ;m
	HR=omega^2*(1+omega^2*tau_rad^2*gam)/g/ang/H0_mod/(1+omega^2*tau_rad^2*gam^2)-1./4/H0_mod^2
	HI=(gam-1)*tau_rad*omega^3/(1+omega^2*tau_rad^2*gam^2)/g/ang/H0_mod

	kr_sq=0.5*(hr+(hr^2+hi^2)^0.5)
	ki_sq=0.5*(-hr+(hr^2+hi^2)^0.5)

	z=2000e3/ang ;m - height of chromosphere
	amp=mdi_pow*exp(z/2/H0)*exp(-1.*sqrt(ki_sq)*z)

    out[*,i,0]=amp

ENDFOR


model=mean(out,dim=2)
cart=plot(mdi_freq*1e3,model[*,1]/max(model[*,1]),color='red',thick=3,xsty=1,ysty=1,yr=[0,1],name='Coronal',$
		  xtitle="Frequency (mHz)",ytitle='Normalised Power',/buffer)
cart2=plot(mdi_freq*1e3,mdi_pow/max(mdi_pow[400:700]),overplot=cart,color='blue',thick=3,name='MDI')
lin=plot([vac,vac]*1e3,[0,1],overplot=cart,thick=3,linestyle=2,name='Acoustic cut-off')
leg=legend(target=[cart,cart2,lin],position=[0.9,0.8],/norm)
cart.save,outpath+'toy_model.png'
stop
END