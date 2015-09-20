@time_lib
@~/dusty/transport/contour_advection/boundary_files

date=convert_date('2002/9/4-12')

datestr1=string(date.year, date.month, date.day, format="(i4.4, i2.2, i2.2)")
if date.hour ne 0 then datestr1=datestr1+"-"+string(date.hour, format="(i2.2)")

datafile="../RET"+datestr1+".0.001sh.36S.hr.1.vec"
outfile="test4-12.ps"

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

spawn, "ls c???.result.bev", contour_files

get_ecmwf_field, date, "SQ", sq, lon, lat, /reverse
nlon=n_elements(lon)
sqfield=[sq[nlon/2:nlon-1, *, 36], sq[0:nlon/2-1, *, 36]]
lon=[lon[nlon/2:nlon-1]-360, lon[0:nlon/2-1]]

set_plot, "z"
device, decomposed=0
;device, filename=outfile

device, set_resolution=[700, 400]

;device, /landscape, xsize=9, ysize=6.5, xoffset=1, yoffset=10, /inches
;device, /color, bits_per_pixel=8

r=indgen(21)*10+40
g=fltarr(21)+40
b=250-indgen(21)*10

tvlct, r, g, b, 1

;goto, legend

grey=200
white=255
black=0;
line=25
title="WV isoline, sigma=36, SH=0.001 for "+time_string(date)
map_set, /cyl, title=title
contour, abs(con), lonret, latret, levels=findgen(21)/20, $
		title=title, xtitle="lon", ytitle="lat", $
		c_color=indgen(21)+1, /cell_fill, xstyle=1, ystyle=1, $
		/overplot, charsize=1.

image1=tvrd()

map_set, /cyl, title=title
for i=0L, n_elements(contour_files)-1 do begin
  tindex=scan_boundary_file(contour_files[i])
  ind=bin_search(tindex, date, /lower)
  get_boundary, contour_files[i], ind, date1, x, y
  ind1=where(x ge 180, cnt)
  if cnt ne 0 then x[ind1]=x[ind1]-360
  if time_compare(date, date1) ne 0 then begin
    Print, time_string(date), time_string(date1)
  endif
  oplot, x, y, thick=1.5, color=line
endfor

;contour, sqfield, lon, lat, levels=0.001, thick=5, /overplot

image2=tvrd()
ind=where(image2 eq line)
image1[ind]=black
tvlct, r, g, b, /get
write_tiff, "test4-12.tif", reverse(image1, 2), red=r, green=g, blue=b

legend:
device, set_resolution=[480, 100]

;device, /close_file
;device, filename="legend.ps"
;device, /portrait, xsize=6.5, ysize=1.5, /inches
x=findgen(21)/40+0.5
contour, rebin(x, 21, 2), x, [0, 1], ystyle=4, $
		title="conditional probability", levels=x, /fill, $
		c_color=indgen(21)+1, charsize=1.

;device, /close_file
;device, /close
image=tvrd()
write_tiff, "legend.tif", reverse(image, 2), red=r, green=g, blue=b

set_plot, "x"

end
