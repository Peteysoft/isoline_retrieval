@rawread_amsu_l1c

;path1="../data/cret/20030101/"
path2="../data/cret/20030102/"

;spawn, "ls "+path1, filelist1, count=nf1
spawn, "ls "+path2, filelist2, count=nf2

set_plot, "ps"
device, /landscape, xsize=9.5, ysize=6.5, xoffset=1, yoffset=10.5, /inches
device, filename="testret4.ps"
device, /color, bits_per_pixel=8

plot, findgen(360)-180, findgen(360)/2-90, /nodata

;ind1=where(strmid(filelist1, 6, 1) eq 'M', cnt1)
;ind2=where(strmid(filelist2, 6, 1) eq 'M', cnt2)

r=findgen(21)*10
g=fltarr(21)
b=250-findgen(21)*10

tvlct, r, g, b, 1

cnt1=0
for i=0L, cnt1-1 do begin
  data=rawread_amsu_l1c(path1+filelist1[ind1[i]])
  ind=where(data.bt gt 0)
  oplot, (data.lon)[ind], (data.lat)[ind], psym=3
endfor

for i=0L, nf2-1 do begin
  data=rawread_amsu_l1c(path2+filelist2[i], missing=missing)
  ind=where(data.bt gt missing)
  ;round_sym, 0.1, /fill
  plots, (data.lon)[ind], (data.lat)[ind], psym=2, color=fix((data.bt)[ind]*10+10)+1, symsize=0.002
  ;break
endfor

date='2003/1/2-12'
get_ecmwf_field, date, "SQ", sq, lon, lat, /reverse
nlon=n_elements(lon)
sqfield=[sq[nlon/2:nlon-1, *, 36], sq[0:nlon/2-1, *, 36]]
lon=[lon[nlon/2:nlon-1]-360, lon[0:nlon/2-1]]

contour, sqfield, lon, lat, level=0.001, thick=3, /overplot

device, /close_file
device, /close
set_plot, "x"

end

