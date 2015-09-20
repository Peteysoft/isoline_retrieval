@boundary_geometry

;set the bins:
nbins=100
minbin=0.00001
maxbin=0.005
dbin=(maxbin/minbin)^(1./(nbins-1))
binvalues=minbin*dbin^findgen(nbins)

goto, skip
res=1./5

nlon=360/res
nlat=180/res+1

lon=findgen(nlon)*res-180
lat=findgen(nlat)*res-90
dlon=res
dlat=res
dlon_ec=1.5
dlat_ec=1.5
nlon_ec=360./dlon_ec
nlat_ec=180./dlat_ec+1

slev=36		;sigma level
SQ_thresh=0.001 ;isoline

;file and date:

nf=10
datelist=time_array(convert_date('2002/9/3'), convert_date('12:0'), nf)
filelist="../RET"+['20020903', $
	'20020903-12', $
	'20020904', $
        '20020904-12', $
        '20020905', $
        '20020905-12', $
        '20020906', $
        '20020906-12', $
        '20020907', $
        '20020907-12']+".0.001sh.36S.hr.1.vec"

;delvar, sq
;delvar, con

station_list=lonarr(10000)
date_list=replicate({time_str}, 10000)

;set the bins:
nbins=30
minbin=0.000001
maxbin=0.01
dbin=(maxbin/minbin)^(1./(nbins-1))
binvalues=minbin*dbin^findgen(nbins)
bins1=lonarr(nbins)
bins2=lonarr(nbins)

maxn=1000000
sqlist=fltarr(maxn)
n=0

for i=0, nf-1 do begin
  print, i
  file=filelist[i]
  date=datelist[i]

  ;get the retrieved conditional probabilities:
  nvar=0L
  r=fltarr(nlon, nlat)
  openr, lun, file, /get_lun
  readu, lun, nvar, r
  free_lun, lun

  ;get the ECMWF field:
  get_ecmwf_field, date, "SQ", sq, eclon, eclat, /reverse
  sqfield=sq[*, *, slev]
  sqfield=[sqfield, sqfield[0, *]]
  eclon=[eclon, 360]

  ;get the retrieved isoline:
  contour, r, lon, lat, levels=0, path_info=path_info, path_xy=path_xy, $
		/path_data_coords

  ;interpolate ecmwf at points along the isoline and for a histogram:
  for j=0, n_elements(path_info)-1 do begin
    l1=path_info[j].offset
    n1=path_info[j].n
    l2=l1+n1-1
    lonpath=path_xy[0, l1:l2]
    latpath=path_xy[1, l1:l2]

    ind=where(lonpath lt 0., cnt)
    if cnt ne 0 then lonpath[ind]=360.+lonpath[ind]

;    ds=point_spacing(lonpath, latpath, /nowrap)
;    s=[0, total(ds, /cum)]

;    if s[n-1] lt 10000. then continue
;    sgrid=[findgen(s[n1-1]/10000.)*10000, s[n1-1]]

;    x=interpol(lonpath, s, sgrid)
;    y=interpol(latpath, s, sgrid)

;    xind=interpol(findgen(nlon_ec+1), eclon, x)
;    yind=interpol(findgen(nlat_ec), eclat, y)
;    sqpath=interpolate(sqfield, xind, yind)
;    binloc=value_locate(binvalues, sqpath)
;    bins1[binloc]=bins1[binloc]+1

    xind=interpol(findgen(nlon_ec+1), eclon, lonpath)
    yind=interpol(findgen(nlat_ec), eclat, latpath)
    sqpath=interpolate(sqfield, xind, yind)
    sqlist[n]=reform(sqpath)
    n=n+n1
    binloc=value_locate(binvalues, sqpath)
    bins2[binloc]=bins2[binloc]+1

;    stop
  endfor
    
  ;lets see...

  ;stop

endfor

sqlist=sqlist[0:n-1]

skip:

binloc=value_locate(binvalues, sqlist)
bins=histogram(binloc, nbins=nbins, min=0, max=nbins)
;bins=histogram(binloc)

xbin=exp((alog(binvalues[1:nbins-1])+alog(binvalues[0:nbins-2]))/2)

ave=mean(sqlist)
std=stddev(sqlist)

xtitle=textoidl("specific humidity [kg kg^{-1}]")

set_plot, "ps"
device, filename="ecmwf_val.ps"
device, /landscape, ysize=6.5, xsize=9, xoffset=1, yoffset=10, /inches

plot_oi, xbin, bins, title="Histogram of ECMWF specific humidity along retrieved isolines",$
		xtitle=xtitle, ytitle="n", thick=3, charsize=1.4
oplot, [sq_thresh, sq_thresh], !y.crange
oplot, [ave, ave], !y.crange, linestyle=2
oplot, [ave-std, ave-std], !y.crange, linestyle=1
oplot, [ave+std, ave+std], !y.crange, linestyle=1

xleg=[2e-5, 3e-5]
oplot, xleg, [35000, 35000]
oplot, xleg, [33000, 33000], linestyle=2
oplot, xleg, [31000, 31000], linestyle=1

xyouts, 4e-5, 35000, "isoline", /data, charsize=1.3
xyouts, 4e-5, 33000, "mean", /data, charsize=1.3
xyouts, 4e-5, 31000, "std. dev.", /data, charsize=1.3

device, /close_file
device, /close
set_plot, "x"

end

