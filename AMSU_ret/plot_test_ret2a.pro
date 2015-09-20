datafile="RET20020904.new.vec"
datafile="dum.vec"
date='2002/9/4'
;outfile="testret20020904.new.ps"
outfile="check.ps"

nlon=360*2
nlat=361

openr, 1, datafile
nvar=0L
con=fltarr(nlon, nlat)
readu, 1, nvar, con
close, 1

lonret=findgen(nlon)/2-180
latret=findgen(nlat)/2-90.

get_ecmwf_field, date, "SQ", sq, lon, lat, /reverse
nlon=n_elements(lon)
sqfield=[sq[nlon/2:nlon-1, *, 36], sq[0:nlon/2-1, *, 36]]
lon=[lon[nlon/2:nlon-1]-360, lon[0:nlon/2-1]]

set_plot, "ps"
device, filename=outfile
device, /landscape, xsize=9, ysize=6.5, xoffset=1, yoffset=10, /inches
device, /color, bits_per_pixel=8

tvlct, [200, 150], [0, 150], [0, 150], 1

grey=2
red=1
white=255
black=0;
title="WV isoline, sigma=36, SH=0.001 for "+date
contour, con, lonret, latret, levels=[-0.85, 0, 0.85], $
		title=title, xtitle="lon", ytitle="lat", $
		c_color=[grey, grey, white], /fill
contour, con, lonret, latret, levels=0, $
		c_thick=5, /overplot

contour, sqfield, lon, lat, levels=0.001, thick=5, /overplot, c_color=ret

device, /close_file
device, /close
set_plot, "x"

end
