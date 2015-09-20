@time_lib
@~/dusty/transport/contour_advection/boundary_files


set_plot, "ps"


resolution=1./5
nlon=360/resolution
nlat=180/resolution+1

lonret=findgen(nlon)*resolution-180
latret=findgen(nlat)*resolution-90.

spawn, "ls c???.result.bev", contour_files

datelist=time_array(convert_date('2002/9/3-12'), convert_date('0/0/1'), 5)

for it=1, 1 do begin

date=datelist[it]

datestr1=string(date.year, date.month, date.day, format="(i4.4, i2.2, i2.2)")
if date.hour ne 0 then datestr1=datestr1+"-"+string(date.hour, format="(i2.2)")

datafile="../RET"+datestr1+".0.001sh.36S.hr.1.vec"
outfile="ret"+datestr1+".col.ps"

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

!p.multi=[0, 1, 2]

;r=[indgen(21)*155/20+100, 70]
;g=r
;b=r

r=[255-indgen(21)*255/20, 150, 255]
g=[fltarr(21), 150, 255]
b=[indgen(21)*255/20, 150, 255]

tvlct, r, g, b, 1
color=21-indgen(21)

;light_grey=11
light_grey=22
red=21
;dark_grey=22
white=23

;goto, legend

cadv_thick=2

title="WV isoline, sigma=36, SH=0.001 for "+time_string(date)

contour, con, lonret, latret, levels=[-0.8, 0.8], c_color=[light_grey, white], /fill, $
	title=title, xtitle="lon", ytitle="lat", xstyle=1, ystyle=1
plot_continents
contour, con, lonret, latret, levels=0, c_thick=2, /overplot, c_color=red
contour, sqfield, lon, lat, levels=0.001, c_thick=3, $
		/overplot

oplot, [50, 65], [-48, -48], thick=2, color=red
oplot, [50, 65], [-58, -58], thick=3
polyfill, [50, 65, 65, 50], [-70, -70, -65, -65], color=light_grey, /data

xyouts, 70, -50, "retrieved", /data
xyouts, 70, -60, "ECMWF", /data
xyouts, 70, -70, "90% confidence"

;!p.multi=[1, 1, 2]
title="Advected contour, sigma=36, SH=0.001 for "+time_string(date)
contour, abs(con), lonret, latret, levels=findgen(21)/20, $
		title=title, xtitle="lon", ytitle="lat", $
		c_color=color, /fill, xstyle=1, ystyle=1

for i=0L, n_elements(contour_files)-1 do begin
  tindex=scan_boundary_file(contour_files[i])
  ind=bin_search(tindex, date, /lower)
  get_boundary, contour_files[i], ind, date1, x, y
  ind1=where(x ge 180, cnt)
  if cnt ne 0 then x[ind1]=x[ind1]-360
  if time_compare(date, date1) ne 0 then begin
    print, time_string(date), time_string(date1)
  endif
  npt=n_elements(x)
  bind=where(abs(x[0:npt-2]-x[1:npt-1]) ge 180, cb)
  if cb gt 0 then begin
    oplot, x[0:bind[0]], y[0:bind[0]], thick=cadv_thick
    for j=1, cb-1 do begin
      oplot, x[bind[j-1]+1:bind[j]], y[bind[j-1]+1:bind[j]], thick=cadv_thick
    endfor
    oplot, x[bind[cb-1]+1:npt-1], y[bind[cb-1]+1:npt-1], thick=cadv_thick
    ;stop
  endif else begin
    oplot, x, y, thick=cadv_thick
  endelse
endfor

;contour, sqfield, lon, lat, levels=0.001, thick=5, /overplot

legend:

device, /close_file
device, /close

;stop

endfor

!p.multi=[0, 1, 1]

device, filename="legend.ps"
device, /portrait, xsize=6.5, ysize=1., /inches
x=findgen(21)/20
contour, rebin(x, 21, 2), x, [0, 1], ystyle=4, $
		title="confidence", levels=x, /fill, $
		c_color=color, xmargin=[3, 3], ymargin=[2, 2]

device, /close_file
device, /close
set_plot, "x"


end
