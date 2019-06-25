
;Loops over solar radii to calculate
;line ratios as a function of height
;
;Requires human interaction

PRO get_den_rat

;From solar radii 1-1.4
h=1.+findgen(40)/10.

FOR i=0,39 DO BEGIN
	density_ratios,'fe_13',10740,10790,7,13,density,ratio,temp=1.6e6, $
		             radtemp=6000, rphot=h[i]

    IF n_elements(rat_arr) EQ 0 THEN rat_arr=ratio ELSE rat_arr=[[rat_arr],[ratio]

ENDFOR

save,density,rat_arr,h,filename='analysis/comp/line_ratio/fe_13_rat.sav'

END