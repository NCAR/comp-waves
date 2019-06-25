pro comp_gf,I1,I2,I3

dp2=dblarr(620,620)
wd2=dblarr(620,620)
ic2=dblarr(620,620)

for t=0,49 do begin
for i=0,619 do begin
for j=0,619 do begin
analytic_gauss_fit_fast,[I1[i,j,t],I2[i,j,t],I3[i,j,t]],33.,p,x,y
dp2[i,j]=p
wd2[i,j]=x
ic2[i,j]=y
end
end
end


END
