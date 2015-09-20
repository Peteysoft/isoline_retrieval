@boundary_geometry
@time_lib

res=5.

nlon=360*res
nlat=180*res+1

lon=findgen(nlon)/res-180
lat=findgen(nlat)/res-90

;do a tedious line-integral to get the confidence interval:
vmr_thresh=0.001
slev=36

vmrt_string=string(vmr_thresh, format='(f5.3)')
slev_string=string(slev, format='(i2.2)')

t0=convert_date('2002/9/3')
half_day=convert_date('12:0')
date=time_array(t0, half_day, 10)

datestr=string(date.day+100L*(date.month+100L*date.year), format="(i8.8)")
ind=where(date.hour ne 0)
datestr[ind]=datestr[ind]+"-"+string(date[ind].hour, format="(i2.2)")

retfile="../RET"+datestr+"."+vmrt_string+"sh.36S.hr.1.vec"

;stop

nfile=n_elements(date)

confield=fltarr(nlon, nlat)

path_con=0
ds=0

nvar=0L

stop

for fi=0L, nfile-1 do begin
  ;read in the data:
  get_ecmwf_field, date[fi], "SQ", sqfield, ecmwf_lon, ecmwf_lat
  vmr=sqfield[*, *, slev]

  print, "Reading file: ", retfile[fi]
  openr, lun, retfile[fi], /get_lun
  readu, lun, nvar, confield
  free_lun, lun
  confield=abs(confield)

  contour, vmr, ecmwf_lon, ecmwf_lat, $
		levels=vmr_thresh, /path_data_coords, path_xy=path_xy, $
	  	path_info=path_info, /data

  npts=n_elements(path_xy)/2

  ind=where(path_xy[0, *] gt 180)
  path_xy[0, ind]=path_xy[0, ind]-360
  lonpath=interpol(findgen(nlon), lon, path_xy[0, *])

  latpath=interpol(findgen(nlat), lat, path_xy[1, *])

  path_con1=interpolate(confield, lonpath, latpath)
  path_con1=reform(path_con1)
  if n_elements(path_con) eq 0 then begin
    path_con=path_con1
  endif else begin
    path_con=[path_con, path_con1]
  endelse

  ninfo=n_elements(path_info)

  ds2=fltarr(npts)
  for i=0L, ninfo-1 do begin
    j1=path_info[i].offset
    n=path_info[i].n
    j2=j1+n-1
    ds1=point_spacing(lonpath[j1:j2], latpath[j1:j2], /nowrap)
    if n gt 2 then begin
      ds2[j1]=[ds1[0]/2, (ds1[0:n-3]+ds1[1:n-2])/2, ds1[n-2]/2]
    endif else begin
      ds2[j1]=[ds1[0]/2, ds1[n-1]/2]
    endelse
  endfor

  if n_elements(ds) eq 0 then begin
    ds=ds2
  endif else begin
    ds=[ds, ds2]
  endelse

;  stop

endfor

nt=n_elements(ds)

s=total(ds)
sind=sort(path_con)
dss=ds[sind]
  
conint=fltarr(nt)
conint[0]=dss[0]
for i=1L, nt-1 do conint[i]=conint[i-1]+dss[i]
conint=conint/s

path_cons=path_con[sind]

set_plot, "ps"
device, filename="conintfunc2.ps"
device, /landscape, xsize=9, ysize=6.5, xoffset=1, yoffset=10, /inches
device, /Helvetica

plot, path_cons, conint, xtitle="confidence rating", $
	  	ytitle="confidence interval", yrange=[0, 1], $
		title="Isoline vs. classification confidence", charsize=1.5, $
		thick=1.5

;interpolate 50%, 90% and 1 std. dev. intervals:
c90=interpol(path_cons, conint, 0.9)
c50=interpol(path_cons, conint, 0.5)
cstd=interpol(path_cons, conint, 0.68269)

oplot, [0, c90], [0.9, 0.9]
oplot, [c90, c90], [0, 0.9]
oplot, [0, c50], [0.5, 0.5]
oplot, [c50, c50], [0, 0.5]
oplot, [0, cstd], [0.68269, 0.68269]
oplot, [cstd, cstd], [0, 0.68260]

xyouts, 0.05, 0.91, '90 %', /data, charsize=1.3
xyouts, 0.05, 0.51, '50 %', /data, charsize=1.3
xyouts, 0.05, 0.69269, '1 std. dev.', /data, charsize=1.3

device, /close_file
device, /close
set_plot, "x"

save, path_cons, conint, filename="confunc2.idlsav"

end

