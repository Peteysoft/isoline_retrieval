@time_lib
@~/dusty/transport/contour_advection/boundary_files


set_plot, "ps"

;goto, legend

resolution=1./5
nlon=360/resolution
nlat=180/resolution+1

lonret=findgen(nlon)*resolution-180
latret=findgen(nlat)*resolution-90.

spawn, "ls c???.result.bev", contour_files

date=convert_date('2002/9/4-12')

datestr1=string(date.year, date.month, date.day, format="(i4.4, i2.2, i2.2)")
if date.hour ne 0 then datestr1=datestr1+"-"+string(date.hour, format="(i2.2)")

datafile="../RET"+datestr1+".0.001sh.36S.hr.1.vec"
outfile="cret_detail.ps"

openr, 1, datafile
nvar=0L
con=fltarr(nlon, nlat)
readu, 1, nvar, con
close, 1

get_ecmwf_field, date, "SQ", sq, lon, lat, /reverse
nlon_ec=n_elements(lon)
sqfield=[sq[nlon_ec/2:nlon_ec-1, *, 36], sq[0:nlon_ec/2-1, *, 36]]
lon=[lon[nlon_ec/2:nlon_ec-1]-360, lon[0:nlon_ec/2-1]]

device, filename=outfile
device, /portrait, xsize=6.5, ysize=6.5, /inches
device, /color, bits_per_pixel=8

r=[indgen(21)*10+55, 100]
g=r
b=r

tvlct, r, g, b, 1
color=21-indgen(21)

light_grey=10
dark_grey=22
white=21

cadv_thick=2

title="WV isoline, sigma=36, SH=0.001 for "+time_string(date)

contour, con, lonret, latret, levels=[-0.7, -0.3], /fill, c_color=[light_grey, dark_grey], $
	xrange=[-95, -65], yrange=[45, 70], $
	title=title, xtitle="lon", ytitle="lat", xstyle=1, ystyle=1, c_label=[1, 1, 1]
contour, con, lonret, latret, levels=[-0.7, -0.3], c_label=[1, 1], /overplot
plot_continents
;contour, con, lonret, latret, levels=0, c_thick=2, /overplot, c_color=dark_grey
;contour, sqfield, lon, lat, levels=0.001, c_thick=3, /overplot

round_sym, 0.3, /fill

for i=0L, n_elements(contour_files)-1 do begin
  tindex=scan_boundary_file(contour_files[i])
  ind=bin_search(tindex, date, /lower)
  get_boundary, contour_files[i], ind, date1, x, y
  ind1=where(x ge 180, cnt)
  if cnt ne 0 then x[ind1]=x[ind1]-360
  if time_compare(date, date1) ne 0 then begin
    print, time_string(date), time_string(date1)
    stop
  endif
  npt=n_elements(x)
  bind=where(abs(x[0:npt-2]-x[1:npt-1]) ge 180, cb)
  if cb gt 0 then begin
    oplot, x[0:bind[0]], y[0:bind[0]], thick=cadv_thick, psym=-8
    for j=1, cb-1 do begin
      oplot, x[bind[j-1]+1:bind[j]], y[bind[j-1]+1:bind[j]], thick=cadv_thick, psym=-8
    endfor
    oplot, x[bind[cb-1]+1:npt-1], y[bind[cb-1]+1:npt-1], thick=cadv_thick, psym=-8
    ;stop
  endif else begin
    oplot, x, y, thick=cadv_thick, psym=-8
  endelse
endfor

oplot, [-93, -92, -91], [67, 67, 67], thick=cadv_thick, psym=-8

xyouts, -90, 67, "advected contour", /data
;contour, sqfield, lon, lat, levels=0.001, thick=5, /overplot

device, /close_file
device, /close

legend:

device, filename="detail_legend.ps"

device, /portrait, xsize=6.5, ysize=1.5, /inches

x=-findgen(11)/10
l=rebin(x, 11, 2)
contour, l, x, [0, 1], levels=[-0.7, -0.3], ystyle=4, title="R", $
		c_color=[light_grey, dark_grey], /fill, $
		xtickv=[-1, -0.7, -0.3, 0], xticks=3

device, /close_file
device, /close

set_plot, "x"


end
