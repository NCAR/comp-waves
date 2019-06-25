;PURPOSE - Clean up comp velocity data! Certain CoMP data sets suffer from
;          aerosol streaks in certain frames and dead pixels
;
;HOW -  Manually provide frame numbers in wave_tracking.pro via ftc keyword
; 		to remove and replace with time median data. 
;       Then uses multi-scale approach to despike
;
;INPUT date - date of observations to clean
;
;OPTIONAL INPUTS ftc - frame numbers that contain streaks
;				 sigma_cut - level at which to remove data spikes
;							 default is 3 sigma
;
;OUPUTS - out cleaned data
;
;HISTORY - Created by RJ Morton 01/2017
;
;
;


PRO comp_clean,date,ftc=ftc,sigma_cut=sigma_cut

inpath='analysis/comp/wave_tracking_output/' 

restore,inpath+date+'/cube_ivw_'+date+'.sav'

;width of median filter
width=7

data=cube_v
sz=size(data)
nx=sz(1)
ny=sz(2)
nt=sz(3)
out=data
temp_ind=index

;removes large scale artefacts on specified frames
;Works by removing poor frames and replacing them
;with time median value of surrounding frames
med=data
nclean=n_elements(ftc)

IF nclean GT 0 THEN BEGIN
	FOR i=0,nx-1 DO FOR j=0,ny-1 DO IF index.mask[i,j] EQ 1 THEN med[i,j,*]=median(reform(data[i,j,*]),width)

	diff=data-med
	mn=moment(diff,dim=3)

    print,'Removing bad frames'
	FOR i=0,nclean-1 DO BEGIN
		out[*,*,ftc[i]]=med[*,*,ftc[i]]

		counter,i,nclean
	ENDFOR
	temp_ind=add_tag(temp_ind,ftc,'Frames_cleaned')


ENDIF

;Removes spikes from images
;Intended to remove dead pixels but also removes large
;noise spikes
;Uses multi-scale filtering to isolate noise following
;Starck & Murtagh
;then determines which pixels are sigma_cut* stand.dev. greater
;than mean value 
;replaces with median of 7 by 7 neighbouring pixels
;
all_spikes=intarr(nx,ny,nt)
;Despike images
IF NOT keyword_set(sigma_cut) THEN sigma_cut=3
print,'Removing spikes'
FOR i=0,nt-1 DO BEGIN
		FOR j=0,1 DO BEGIN
	        image2=out[0:-1,0:-1,i]
	        spikeim=intarr(nx,ny)

	        atrous,image2,decomp=wav,n_scales=1
	        wav1=wav[0:-1,0:-1,1]
	        
	        ;find distribution for values where mask is one
	        mn=moment(wav1[where(index.mask eq 1)])

	        smed=median(wav1,7)
            in=where(wav1 gt sigma_cut*sqrt(mn[1]) or wav1 lt -sigma_cut*sqrt(mn[1]))
	        spikeim[in]=1
	        all_spikes[0:-1,0:-1,i]=all_spikes[0:-1,0:-1,i]+spikeim

	        ;replace data and recombine
	        wav1[in]=smed[in]
	  		out[0:-1,0:-1,i] = wav[0:-1,0:-1,0]+wav1
   		ENDFOR
   	 counter,i,nclean	
ENDFOR
index_out=add_tag(temp_ind,'yes','Despiked')

cube_v=out
index=index_out

save,cube_i,cube_v,cube_w,index,filename=inpath+date+'/cube_ivw_'+date+'.sav'


;FOR k=0,1 do BEGIN
;IF k eq 0 THEN temp=data ELSE temp=out
;IF k eq 0 THEN sd=3. ELSE sd=5.
;count=!null

;FOR i=0,nx-1 DO FOR j=0,ny-1 DO IF index.mask[i,j] EQ 1 THEN med[i,j,*]=median(reform(temp[i,j,*]),width)

;diff=data-med
;mn=moment(diff,dim=3)


;FOR i=0,610 DO FOR j=0,619 DO BEGIN
;   	IF index.mask[i,j] EQ 1 THEN BEGIN
;		samp=reform(diff[i,j,*])
;		in=where(samp gt sd*sqrt(mn[i,j,1]))
;		in2=where(samp lt -sd*sqrt(mn[i,j,1]))
;	    IF in[0] NE -1 THEN IF in2[0] NE -1 THEN BEGIN 
;	       out[i,j,[in,in2]]=med[i,j,[in,in2]]
;	       IF n_elements(count) EQ 0 THEN count=n_elements(in)+n_elements(in2) ELSE $
;	       count=[count,n_elements(in)+n_elements(in2)]
;	    ENDIF   
;	    IF in[0] NE -1 THEN IF in2[0] EQ -1 THEN BEGIN
;	    	out[i,j,[in]]=med[i,j,[in]]
;	    	IF n_elements(count) EQ 0 THEN count=n_elements(in) ELSE $
;	       count=[count,n_elements(in)]
;	    ENDIF
;	    IF in[0] EQ -1 THEN IF in2[0] NE -1 THEN BEGIN
;	    	out[i,j,[in2]]=med[i,j,[in2]]
;	    	IF n_elements(count) EQ 0 THEN count=n_elements(in2) ELSE $
;	       count=[count,n_elements(in2)]
;	    ENDIF
;   ENDIF
;ENDFOR
;ENDFOR

END