@boundary_geometry

nlon=240L
nlat=121L

;do a tedious line-integral to get the confidence interval:
vmr_thresh=0.001
slev=36

vmrt_string=string(vmr_thresh, format='(f5.3)')
slev_string=string(slev, format='(i2.2)')

dates=['1998082718', '1999030218', '1999092212', '2000101600', $
		'2001122006', '2002061212']
confile="ecmwf"+dates+'_cld'+slev_string+'s'+vmrt_string+"SQ.ret.b.con"
vmrfile="ecmwf"+dates+"SQ"+slev_string+"s.dat"

;stop

nfile=n_elements(dates)

vmr=fltarr(nlon, nlat)
confield=fltarr(nlon, nlat)

path_con=0
ds=0

for fi=0L, nfile-1 do begin
  ;read in the data:
  print, "Reading file: ", vmrfile[fi]
  openr, lun, vmrfile[fi], /get_lun
  readu, lun, vmr
  free_lun, lun

  print, "Reading file: ", confile[fi]
  openr, lun, confile[fi], /get_lun
  readu, lun, confield
  free_lun, lun

  contour, vmr, levels=vmr_thresh, /path_data_coords, path_xy=path_xy, $
	  	path_info=path_info

  npts=n_elements(path_xy)/2

  path_con1=interpolate(confield, path_xy[0, *], path_xy[1, *])
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
    ds1=point_spacing(path_xy[0, j1:j2], path_xy[1, j1:j2], /nowrap)
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
device, filename="conintfunc.ps"
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

end

