@sonde_file_functions

;goto, skip

nlon=360*5
nlat=180*5+1

lon=findgen(nlon)/5-180
lat=findgen(nlat)/5-90
dlon=0.2
dlat=0.2
dlon_ec=1.5
dlat_ec=1.5

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

openw, 1, "station_list2.txt"

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

  ;get the ECMWF pressure field:
  get_ecmwf_field, date, "PP", pp, eclon, eclat

  ind=where(time_compare(index.date, date) eq 0, nsonde)
  printf, 1, time_string(date)
  printf, 1, index[ind].id
  station_list[nstat:nstat+nsonde-1]=index[ind].id
  date_list[nstat:nstat+nsonde-1]=date
  nstat=nstat+nsonde

  ;lets see...

  ;first, interpolate the pressure for each sonde location
  rslon=index[ind].lon
  sub=where(rslon lt 0, cnt)
  if cnt gt 0 then rslon[sub]=360+rslon[sub]
  xind_ec=rslon/dlon_ec
  yind_ec=(index[ind].lat+90)/dlat_ec

  p=interpolate(pp[*, *, slev], xind_ec, yind_ec)*100

  ;interpolate the humidity values for each sonde
  sq1=fltarr(nsonde)
  for j=0, nsonde-1 do begin
    read_sonde_profile, sonde_file, ind[j], pprof, tprof, zprof, qprof
    sq1[j]=interpol(qprof, alog(pprof), alog(p[j]))
  endfor

  ;convert the humidity values and add to list:
  ;if n_elements(sq) eq 0 then sq=0.622*sq1 else sq=[sq, 0.622*sq1]
  sq[nstat:nstat+nsonde-1]=0.622*sq1

  ;interpolate the cond. prob. for each sonde location
  xind=(index[ind].lon+180)/dlon
  yind=(index[ind].lat+90)/dlat
  con1=interpolate(r, xind, yind)
  ;if n_elements(con) eq 0 then con=(1+con1)/2 else con=[con, (1+con1)/2]
  con[nstat:nstat+nsonde-1]=(1+con1)/2.

  ;stop

endfor

close, 1

sq=sq[0:nstat-1]
con=con[0:nstat-1]

skip:

;sort the cond. prob.:
sind=sort(con)

;should be a straight line:
set_plot, "ps"
device, filename="rs_val.ps"
device, /landscape, ysize=4, xsize=9, xoffset=1, yoffset=10, /inches

!p.multi=[0, 2, 1]

cons=con[reverse(sind)]
sqs=sq[reverse(sind)]

ind1=where(cons lt 0.5, cnt1)
y1=total((sqs[ind1] lt sq_thresh)/(1-cons[ind1]), /cum)
x1=findgen(cnt1)

acc1=total((sqs[ind1] gt sq_thresh) eq (cons[ind1] gt 0.5))/cnt1
ytitle1=textoidl("\Sigma_{i=1}^n \delta_{1 c(i)}/P(1|x)")

round_sym, 0.3, /fill

plot, x1, y1, xtitle="n", ytitle=ytitle1, xstyle=1, psym=8, $
		thick=5, subtitle="accuracy ="+string(acc1, format="(f5.3)")
;stop
oplot, findgen(cnt1), findgen(cnt1)

;xtickname=string(findgen(6)/5*0.5+0.5, format="(f3.1)")
;xtickv=interpol(x1, 1-cons[ind1], findgen(6)/5*0.5+0.5)
xtickv=[0, 200, 400, 600, 800, x1[cnt1-1]]
xtickname=string(interpol(1-cons[ind1], x1, xtickv), $
		format="(f4.2)")
axis, 0, 1, /normal, xtickname=xtickname, xtickv=xtickv, $
		xstyle=1, xticks=5
xyouts, 0.25, 0.88, "P(1|x)", /normal

cons=con[sind]
sqs=sq[sind]

ind2=where(cons gt 0.5, cnt2)
y2=total((sqs[ind2] gt sq_thresh)/cons[ind2], /cum)
x2=findgen(cnt2)

acc2=total((sqs[ind2] gt sq_thresh) eq (cons[ind2] gt 0.5))/cnt2
ytitle2=textoidl("\Sigma_{i=1}^n \delta_{2 c(i)}/P(2|x)")

plot, x2, y2, xtitle="n", ytitle=ytitle2, xstyle=1, psym=8, $
		thick=5, subtitle="accuracy ="+string(acc2, format="(f5.3)")

oplot, findgen(cnt2), findgen(cnt2)

xtickv=[0, 20, 40, 60, 100, 120, 140, x2[cnt2-1]]
xtickname=string(interpol(cons[ind2], x2, xtickv), $
		format="(f4.2)")
axis, 0, 1, /normal, xtickname=xtickname, xtickv=xtickv, xstyle=1
xyouts, 0.75, 0.88, "P(2|x)", /normal

device, /close_file
device, /close
set_plot, "x"

!p.multi=[0, 1, 1]

;(this is ineffective--let's do it in the usual manner--by binning):

nbin=5
bin=fix(con*nbin)

conret=fltarr(nbin)
freq=fltarr(nbin)
for i=0, nbin-1 do begin
  ind=where(bin eq i, cnt)
  if cnt eq 0 then begin
    conret[i]=1.0*i/nbin+0.5
    freq[i]=0./0.
    continue
  endif
  conret[i]=total(con[ind])/cnt
  freq[i]=total(sq[ind] gt SQ_thresh)/cnt
endfor

;generate the station list:
sind2=sort(station_list[0:nstat-1])
station_list2=station_list[sind2]
date_list2=date_list[sind2]
uind=uniq(station_list2)

sflag=intarr(nf, n_elements(uind))

for i=0, nf-1 do begin
  ind=where(time_compare(date_list2, datelist[i]) eq 0, cnt)
  for j=0, cnt-1 do begin
    sflag[i, bin_search(uind, ind[j], /upper)]=1
  endfor
endfor

openw, 1, "station_list.txt"

printf, 1, "station id", datelist.day, format="(a10, 10i3)"
for i=0, n_elements(uind)-1 do begin
  rest=""
  for j=0, nf-1 do begin
    if sflag[j, i] eq 1 then rest=rest+"& X" else rest=rest+"&  "
  endfor
  printf, 1, station_list2[uind[i]], rest, "\\"
endfor

close, 1


end

