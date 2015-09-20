@boundary_geometry
@time_lib

res=1.5

nlon=360/res
nlat=180/res+1

lon=findgen(nlon)*res-180
lat=findgen(nlat)*res-90

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

retfile="../RET"+datestr+"."+vmrt_string+"sh.36S.ecmwf.vec"

;stop

nfile=n_elements(date)

con=fltarr(nlon, nlat, nfile)
confield=fltarr(nlon, nlat)
vmr=fltarr(nlon, nlat, nfile)

path_con=0
ds=0

nvar=0L

;stop

nbin=20
r=fltarr(nbin)
hist=lonarr(nbin)
nperbin=lonarr(nbin)

for fi=0L, nfile-1 do begin
  ;read in the data:
  get_ecmwf_field, date[fi], "SQ", sqfield, ecmwf_lon, ecmwf_lat
  vmr[*, *, fi]=sqfield[*, *, slev]

  print, "Reading file: ", retfile[fi]
  openr, lun, retfile[fi], /get_lun
  readu, lun, nvar, confield
  free_lun, lun

  con[*, *, fi]=confield



;  stop

endfor

for j=0, nbin-1 do begin
  ind=where(fix((con/2+0.5)*nbin) eq j, cnt)
  r[j]=r[j]+total(con[ind])
  hist[j]=hist[j]+total(vmr[ind] gt vmr_thresh)
  nperbin[j]=nperbin[j]+cnt
endfor

r=r/nperbin

set_plot, "ps"
device, filename="accuracy_ecmwf.ps"
device, /portrait, xsize=6.5, ysize=5.8, /inches
device, /Helvetica

round_sym, 1;, /fill
plot, r, 2.*hist/nperbin-1, xtitle="R", $
	  	ytitle="frequency (transformed)", $
		title="ECMWF vs. isoline retrieval", charsize=1.5, $
		psym=-8

oplot, findgen(21)/10-1, findgen(21)/10-1, linestyle=2
oplot, [-1, 1], [0, 0]
oplot, [0, 0], [-1, 1]


device, /close_file
device, /close
set_plot, "x"


end

