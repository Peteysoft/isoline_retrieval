@boundary_geometry
@boundary_files
@time_lib

restore, "confunc.idlsav"

path_cons_sim=path_cons
conint_sim=conint

restore, "confunc2.idlsav"

path_cons_ecmwf=path_cons
conint_ecmwf=conint

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

nt0=5
t0=time_array(convert_date('2002/9/3'), convert_date('0/0/1'), nt0)

half_day=convert_date('12:0')

set_plot, "ps"
device, filename="conintfunc3.ps"
device, /portrait, xsize=6.5, ysize=5.8, /inches
device, /Helvetica

plot, path_cons_ecmwf, conint_ecmwf, xtitle="confidence rating", xrange=[0, 1], $
	  	ytitle="path fraction", yrange=[0, 1], $
		title="Path fraction vs. confidence rating", charsize=1.3

oplot, path_cons_sim, conint_sim, thick=10, color=180

for i=0L, nt0-1 do begin

date=time_array(t0[i], half_day, 2)

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

;stop

spawn, "ls c???.result.bev", cfile

for fi=0L, nfile-1 do begin
  ;read in the data:
;  get_ecmwf_field, date[fi], "SQ", sqfield, ecmwf_lon, ecmwf_lat
;  vmr=sqfield[*, *, slev]

  print, "Reading file: ", retfile[fi]
  openr, lun, retfile[fi], /get_lun
  readu, lun, nvar, confield
  free_lun, lun
  confield=abs(confield)

  for j=0, n_elements(cfile)-1 do begin
    tindex=scan_boundary_file(cfile[j])
    tind=bin_search(tindex, date[fi], /lower)
    get_boundary, cfile[j], tind, date1, x, y

    npts=n_elements(x)

    ind=where(x gt 180, cnt)
    if cnt gt 0 then x[ind]=x[ind]-360
    lonpath=interpol(findgen(nlon), lon, x)

    latpath=interpol(findgen(nlat), lat, y)

    path_con1=interpolate(confield, lonpath, latpath)
    path_con1=reform(path_con1)
    if n_elements(path_con) eq 0 then begin
      path_con=path_con1
    endif else begin
      path_con=[path_con, path_con1]
    endelse

    ds1=point_spacing(lonpath, latpath)
    ds2=[ds1[0]+ds1[npts-1], ds1[0:npts-2]+ds1[1:npts-1]]/2

    if n_elements(ds) eq 0 then begin
      ds=ds2
    endif else begin
      ds=[ds, ds2]
    endelse

  endfor
;  stop

endfor

nt=n_elements(ds)

s=total(ds)
sind=sort(path_con)
dss=ds[sind]
  
conint=fltarr(nt)
conint[0]=dss[0]
for j=1L, nt-1 do conint[j]=conint[j-1]+dss[j]
conint=conint/s

path_cons=path_con[sind]

oplot, path_cons, conint, linestyle=1

sub=n_elements(conint)/2
xyouts, path_cons[sub], conint[sub], string(i+1, format="(i2)"), /data

endfor

oplot, [0.05, 0.1], [0.925, 0.925], thick=10, color=180
oplot, [0.05, 0.1], [0.875, 0.875]
oplot, [0.05, 0.1], [0.825, 0.825], linestyle=1

xyouts, 0.15, 0.925, "simulated", /data
xyouts, 0.15, 0.875, "ECMWF", /data
xyouts, 0.15, 0.825, 'advected contours:', /data
xyouts, 0.1, 0.775, '1: 2002/09/03', /data
xyouts, 0.1, 0.725, '2: 2002/09/04', /data
xyouts, 0.1, 0.675, '3: 2002/09/05', /data
xyouts, 0.1, 0.625, '4: 2002/09/06', /data
xyouts, 0.1, 0.575, '5: 2002/09/07', /data

device, /close_file
device, /close
set_plot, "x"

end

