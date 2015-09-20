datafile="RET20020904-12.0.001sh.36S.hr.1.vec"
date='2002/9/4-12'
outfile="ret.colour.ps"

nlon=360*5
nlat=180*5+1

openr, 1, datafile
nvar=0L
con=fltarr(nlon, nlat)
readu, 1, nvar, con
close, 1

lonret=findgen(nlon)/5-180
latret=findgen(nlat)/5-90.

get_ecmwf_field, date, "SQ", sq, lon, lat, /reverse
nlon=n_elements(lon)
sqfield=[sq[nlon/2:nlon-1, *, 36], sq[0:nlon/2-1, *, 36]]
lon=[lon[nlon/2:nlon-1]-360, lon[0:nlon/2-1]]

;set_plot, "ps"
set_plot, "z"
;device, filename=outfile
;device, /landscape, xsize=9, ysize=6.5, xoffset=1, yoffset=10, /inches
;device, /color, bits_per_pixel=8
device, decomposed=0

device, set_resolution=[700, 400]

tvlct, [0, 255, 200, 170], [0, 255, 0, 170], [0, 255, 0, 170]

red=2
grey=3
white=1
black=0;
title="WV isoline, sigma=36, SH=0.001 for "+date
retthick=2
ecmwfthick=2
contour, con, lonret, latret, levels=[-0.85, 0, 0.85], $
		title=title, xtitle="lon", ytitle="lat", $
		c_color=[grey, grey, white], /fill, $
		xrange=[-180, 180], yrange=[-90, 90], xstyle=1, ystyle=1, $
		charsize=1., background=white, color=black

plot_continents, color=black
contour, sqfield, lon, lat, levels=0.001, thick=ecmwfthick, /overplot, c_color=red

contour, con, lonret, latret, levels=0, color=black, $
		c_thick=retthick, /overplot

oplot, [70, 90], [-50, -50], thick=retthick, color=black
oplot, [70, 90], [-60, -60], thick=ecmwfthick, color=red
polyfill, [70, 90, 90, 70], [-72, -72, -68, -68], color=grey

xyouts, 100, -51, "retrieved", charsize=1., color=black
xyouts, 100, -61, "ECMWF", charsize=1., color=black
xyouts, 100, -71, "90% confidence", charsize=1., color=black

image=tvrd()

tvlct, r, g, b, /get

write_tiff, "ret.colour.tif", reverse(image, 2), red=r, green=g, blue=b

;device, /close_file
;device, /close
set_plot, "x"

end
