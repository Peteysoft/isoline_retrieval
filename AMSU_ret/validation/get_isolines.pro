;gets the isolines of a field and writes them to a bunch of files
;for use in a contour advection simulation

magic=-55L

file="../RET20020903.0.001sh.hr.vec"
date=convert_date('2002/9/3')

res=5.

dlon=1./res
dlat=1./res
nlon=360/dlon
nlat=180/dlat+1

lon=findgen(nlon)*dlon-180
lat=findgen(nlat)*dlat-90

field=fltarr(nlon, nlat)
nvar=1L

openr, 1, file
readu, 1, nvar, field
close, 1

stop

;so the contour advection codes don't choke, smooth and resample the
;isolines at a lower resolution:
;c=linfilt_coeff_gauss(2)
;field=linfilt_2d(field, c)

f2=[field, field]
lon2=[lon, lon+360]

;get the contours:
contour, f2, lon2, lat, levels=0, /path_data_coord, path_xy=path_xy, path_info=path_info

window, 0
contour, field, lon, lat, levels=0

window, 1
plot, findgen(10), xrange=[-180, 180], yrange=[-90, 90], /nodata

;now go through and make sure:
;1. none of the contours touch the end of the map
;2. there are no duplicates

ind=where(path_info.n ge 50, ncont)

done=bytarr(ncont)

for i=0, ncont-1 do begin
  offset=path_info[ind[i]].offset
  n=path_info[ind[i]].n

  dum=where(path_xy[0, offset:offset+n-1] le -180., cnt)
  if cnt gt 0 then continue
  dum=where(path_xy[0, offset:offset+n-1] ge 720.-180.5, cnt)
  if cnt gt 0 then continue

  cflag=0
  clat=path_xy[1, offset:offset+n-1]
  for j=i-1, 0, -1 do begin
    clat2=path_xy[1, path_info[ind[j]].offset:path_info[ind[j]].offset+path_info[ind[j]].n]
    if max(clat) eq max(clat2) and min(clat) eq min(clat2) and done[j] eq 1 then begin
      cflag=1
      break
    endif
  endfor
  if cflag then continue

  clon=(path_xy[0, offset:offset+n-1]+180) mod 360 -180
  oplot, clon, clat, psym=3, color=rgbtoindex(200, 0, 0)
  ;stop
  ind1=where(clon lt 0, cnt)
  if cnt gt 0 then clon[ind1]=360+clon[ind1]
  openw, 1, "c"+string(i, format="(i3.3)")+".bev"
  writeu, 1, magic, 1L, date, n, clon, clat
  close, 1

  done[i]=1

;  plot, clon, clat

;  stop

endfor

end

