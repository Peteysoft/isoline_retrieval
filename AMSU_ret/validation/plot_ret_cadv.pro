@time_lib
@~/dusty/transport/contour_advection/boundary_files

date=convert_date('2002/9/4-12')

datestr1=string(date.year, date.month, date.day, format="(i4.4, i2.2, i2.2)")
if date.hour ne 0 then datestr1=datestr1+"-"+string(date.hour, format="(i2.2)")

datafile="../RET"+datestr1+".0.001sh.36S.hr.1.vec"
outfile="test2.ps"

resolution=1./5
nlon=360/resolution
nlat=180/resolution+1

openr, 1, datafile
nvar=0L
con=fltarr(nlon, nlat)
readu, 1, nvar, con
close, 1

lonret=findgen(nlon)*resolution-180
latret=findgen(nlat)*resolution-90.

spawn, "ls *.result.bev", contour_files

get_ecmwf_field, date, "SQ", sq, lon, lat, /reverse
nlon=n_elements(lon)
sqfield=[sq[nlon/2:nlon-1, *, 36], sq[0:nlon/2-1, *, 36]]
lon=[lon[nlon/2:nlon-1]-360, lon[0:nlon/2-1]]

set_plot, "ps"
device, filename=outfile
device, /landscape, xsize=9, ysize=6.5, xoffset=1, yoffset=10, /inches

grey=200
white=255
black=0;
title="WV isoline, sigma=36, SH=0.001 for "+time_string(date)
contour, con, lonret, latret, levels=[-0.85, 0, 0.85], $
		title=title, xtitle="lon", ytitle="lat", $
		c_color=[grey, grey, white], /fill, charsize=1.5
contour, con, lonret, latret, levels=0, $
		c_thick=5, /overplot

for i=0L, n_elements(contour_files)-1 do begin
  tindex=scan_boundary_file(contour_files[i])
  ind=bin_search(tindex, date, /lower)
  get_boundary, contour_files[i], ind, date1, x, y
  ind1=where(x ge 180, cnt)
  if cnt ne 0 then x[ind1]=x[ind1]-360
  if time_compare(date, date1) ne 0 then begin
    print, time_string(date), time_string(date1)
  endif
  oplot, x, y, psym=3
endfor

;contour, sqfield, lon, lat, levels=0.001, thick=5, /overplot

device, /close_file
device, /close
set_plot, "x"

end
