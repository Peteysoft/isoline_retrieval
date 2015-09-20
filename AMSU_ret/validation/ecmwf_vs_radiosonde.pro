@sonde_file_functions

;goto, skip

dlon=1.5
dlat=1.5

slev=36		;sigma level
SQ_thresh=0.001 ;isoline

;get the sonde index:
;sonde_file="/rinax/storage/users/pmills/polarstern/polarstern.1.rsd"
sonde_file="/rinax/storage/users/pmills/TEMP_database.rsd"
get_sonde_index, sonde_file, index

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

maxn=10000

station_list=lonarr(maxn)
date_list=replicate({time_str}, maxn)

sq=fltarr(maxn)
sq_ec=fltarr(maxn)

nstat=0

for i=0, nf-1 do begin
  print, i
  file=filelist[i]
  date=datelist[i]

  ;get the ECMWF pressure & humidity fields:
  get_ecmwf_field, date, "PP", pp, lon, lat
  get_ecmwf_field, date, "SQ", sqfield, lon, lat

  ind=where(time_compare(index.date, date) eq 0, nsonde)
  station_list[nstat:nstat+nsonde-1]=index[ind].id
  date_list[nstat:nstat+nsonde-1]=date

  ;lets see...

  ;first, interpolate the pressure for each sonde location
  rslon=index[ind].lon
  sub=where(rslon lt 0, cnt)
  if cnt gt 0 then rslon[sub]=360+rslon[sub]
  xind=rslon/dlon
  yind=(90-index[ind].lat)/dlat

  p=interpolate(pp[*, *, slev], xind, yind)*100

  ;interpolate the humidity values for each sonde
  sq1=fltarr(nsonde)
  for j=0, nsonde-1 do begin
    read_sonde_profile, sonde_file, ind[j], pprof, tprof, zprof, qprof
    sq1[j]=interpol(qprof, alog(pprof), alog(p[j]))
  endfor

  ;convert the humidity values and add to list:
  sq[nstat:nstat+nsonde-1]=0.622*sq1
  ;if n_elements(sq) eq 0 then sq=0.622*sq1 else sq=[sq, 0.622*sq1]

  ;interpolate ecmwf humidity values and add to list:
  sq_ec1=interpolate(sqfield[*, *, slev], xind, yind)
  sq_ec[nstat:nstat+nsonde-1]=sq_ec1
  ;if n_elements(sq_ec) eq 0 then sq_ec=sq_ec1 else sq_ec=[sq_ec, sq_ec1]

  ;stop

  nstat=nstat+nsonde
endfor

sq=sq[0:nstat-1]
sq_ec=sq_ec[0:nstat-1]

plot, sq, sq_ec, psym=2

print, total((sq gt sq_thresh) eq (sq_ec gt sq_thresh))/n_elements(sq)


end

