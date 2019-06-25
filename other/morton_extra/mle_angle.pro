;Maximum likelihood estimate of wave angle
;assuming points are spatially distributed as
;multivariate normal

FUNCTION mle_angle,xvec,yvec,len=len

nel=n_elements(xvec)
mean_cent=[mean(xvec),mean(yvec)]
devs=[[xvec-mean_cent[0]],[yvec-mean_cent[1]]]

;Bias corrected covariance matrix
covar_est=transpose(devs)#(devs)/(nel-1.)

;Compute eigenvectors of covariance matrix
egval=eigenql(covar_est,eigenvectors=egvec)

angle=atan(egvec[1,0]/egvec[0,0])*180./!pi

;Calculate ellipse half lengths
;confidence regions should be len*sqrt(chi^2_{p,alpha})
len=sqrt(egval)

return,angle
END