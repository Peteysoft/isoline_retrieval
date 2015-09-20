@wv_formulas
@time_lib
@/home/home02/pmills/climate/block2/animate
@/home/home02/pmills/transport/contour_advection/boundary_geometry

pro ecmwf_iso_settings, cfile, var, nt
  ;set the parameters for doing the retrieval:
  common ecmwf_iso_parm, arts_path, basename, s1, s2, t0lev, reverse_flag, class_file, variables, ntrain, se
  	;basename	base name for all intermediate files
  	;s1:		upper level for calculating tlapse
	;s2:		lower level for calculating tlapse
	;t0lev		level for base temperature
	;reverse_flag	1 if levels in "class file" are reversed
	;class_file	binary file holding profiles plus integrated bt's
	;variables	the variables to use for the classification--as index from 0-6
	;ntrain		number of training vectors to use--same as in class file if undef.

  s1=33
  s2=43
  t0lev=35
  reverse_flag=1
  se=0.95

  if n_elements(cfile) ne 1 then class_file="data/amsub_ecmwf_sim0.95se.dat" $
  		else class_file=cfile
  if n_elements(var) lt 1 or n_elements(var) gt 6 then begin
    variables=[0, 1, 4, 5, 6]
  endif else begin
    variables=var
  endelse

  if n_elements(nt) eq 1 then ntrain=nt

  basename="ecmwf_iso_data/ecmwf_iso_"
  arts_path="arts_files/"
end

pro ecmwf_field_amsu_bt, date, ngrid=ngrid
  common ecmwf_iso_parm, arts_path, base_name, s1, s2, t0lev, reverse_flag, class_file, variables, ntrain, se

  print, "Simulating the brightness temperatures:"
  print, "Writing the profiles to a file..."

  time=convert_date(date)
  datestr=string(time.year, time.month, time.day, time.hour, $
  		format="(i4.4, i2.2, i2.2, i2.2)")

  rsfile=arts_path+"ecmwf_field"+datestr+".aa"
  ffile=arts_path+"f_mono.aa"
  artsfile=arts_path+"ecmwf_amsu"+datestr+".arts"
  ofile=arts_path+"ecmwf_amsu_bt"+datestr+".aa"

  ;convert field to arts ascii:
  ecmwf_field_to_arts, date, rsfile

  ;get the frequency grid:
  f_mono=amsu_b_frequency_grid(ngrids=ngrid)

  ;write out the frequency grid:
  write_arts_ascii, ffile, reform(f_mono, 1, n_elements(f_mono))

  ;write the control file:
  rs_write_arts_control, artsfile, rsfile, ffile, ofile, emiss=se

  ;run arts:
  print, "Starting arts..."
  spawn, "arts "+artsfile

  ;read the output:
  print, "Reading and integrating the output..."
  bt=read_arts_ascii(ofile)

  ;integrate and write to a more convenient format:
  chan18=total((*bt[0])[*, 0:ngrid*2-1], 2)/ngrid/2
  chan19=total((*bt[0])[*, 2*ngrid:4*ngrid-1], 2)/ngrid/2
  chan20=total((*bt[0])[*, 4*ngrid:6*ngrid-1], 2)/ngrid/2
  chan21=total((*bt[0])[*, 6*ngrid:8*ngrid-1], 2)/ngrid/2
  chan22=total((*bt[0])[*, 8*ngrid:10*ngrid-1], 2)/ngrid/2

  openw, lun, base_name+"chan18_bt"+datestr, /get_lun
  writeu, lun, chan18
  free_lun, lun
  openw, lun, base_name+"chan19_bt"+datestr, /get_lun
  writeu, lun, chan19
  free_lun, lun
  openw, lun, base_name+"chan20_bt"+datestr, /get_lun
  writeu, lun, chan20
  free_lun, lun
  openw, lun, base_name+"chan21_bt"+datestr, /get_lun
  writeu, lun, chan21
  free_lun, lun
  openw, lun, base_name+"chan22_bt"+datestr, /get_lun
  writeu, lun, chan22
  free_lun, lun

end

pro ecmwf_training_matrix, date
  common ecmwf_iso_parm, arts_path, base_name, s1, s2, t0lev, reverse_flag, class_file, variables, ntrain

  print, "Generating the training matrix:"
  print, "Reading the data..."

  if n_elements(filename) ne 1 then ecmwf_iso_settings

  time=convert_date(date)
  datestr=string(time.year, time.month, time.day, time.hour, $
  		format="(i4.4, i2.2, i2.2, i2.2)")

  ;generate the training matrix:
  if n_elements(ntrain) ne 0 then nsamp=ntrain
  read_class_file, class_file, data, nrandom=nsamp

  ;delete those records having the same date as the desired retrieval:
  ind=where(time_compare(data.date, time) ne 0, nsamp)
  data=data[ind]

  ;get the mean temperature profile for the training data:
  print, "Calculating the mean temperature profile..."

  if keyword_set(reverse_flag) then begin
    t0leva=59-t0lev
    s1a=59-s2
    s2a=59-s1
  endif else begin
    t0leva=t0lev
    s1a=s1
    s2a=s2
  endelse

  min_corr=1.
  tlapse=fltarr(nsamp)
  t0=fltarr(nsamp)
  for i=0L, nsamp-1 do begin
    tprof=data[i].profile[s1a:s2a].t
    pprof=data[i].profile[s1a:s2a].p
    tlapse[i]=regress(alog(tprof), alog(pprof))
    t0[i]=total(tprof)/n_elements(tprof)
    tlapse[i]=tlapse[i]/9.8
    ;t0[i]=tprof[5]
    corr=correlate(alog(tprof), alog(pprof))
    if corr lt min_corr then min_corr=corr
  endfor
  print, "Minimum correlation between p and t:", min_corr

  ;t0=data.profile[t0leva].t

  print, "Writing the data to a file..."

  ;find the volume mixing ratio:
  vmr_all=magnus_formula(data.profile.dp)/data.profile.p

  ;also write the volume mixing ratio to a file:
  for si=30, 45 do begin
    vmr=vmr_all[59-si, ind]
    sstr=string(si, format="(i2.2)")+"s"
    openw, lun, base_name+"vmr_train"+sstr+datestr+".dat", /get_lun
    writeu, lun, vmr
    free_lun, lun
  endfor

  ;concatinate the training matrix:
  trainmatrix=transpose([	[t0], $
			[tlapse], $
			[data.btemp[0]], $
			[data.btemp[1]], $
			[data.btemp[2]], $
			[data.btemp[3]], $
			[data.btemp[4]]])

  trainmatrix=trainmatrix[variables, *]

  ;write it to a file:
  openw, lun, base_name+"train"+datestr+".vec", /get_lun
  writeu, lun, n_elements(variables)
  writeu, lun, trainmatrix
  free_lun, lun

;  stop

end

pro ecmwf_iso_retrieve, date, sigma, vmr_thresh, field, confidence, wopt=wopt, knn=knn
  common ecmwf_iso_parm, arts_path, base_name, s1, s2, t0lev, reverse_flag, class_file, variables

  if n_elements(knn) ne 1 then knn=10000
  if n_elements(wopt) ne 1 then wopt=50
;  if n_elements(sigma) ne 1 then sigma=33

  if n_elements(s1) eq 0 then ecmwf_iso_settings

  print, "Retrieving an isoline from simulated brightness temperatures:"

  time=convert_date(date)
  datestr=string(time.year, time.month, time.day, time.hour, $
  		format="(i4.4, i2.2, i2.2, i2.2)")

  sstr=string(sigma, format="(i2.2)")+"s"

  ;generate the file names:
  retbase=base_name+"ret"+datestr
  trainbase=base_name+"train"+datestr

  get_ecmwf_field, date, "TT", tt, lon, lat, /noread
  nlon=n_elements(lon)
  nlat=n_elements(lat)


  if file_test(retbase+".vec", /read, /regular) ne 1 then begin

    print, "Generating the retrieval matrix..."

    ;retrieve the temperature profile:
    get_ecmwf_field, date, "TT", tt, lon, lat
    get_ecmwf_field, date, "PP", pp, lon, lat
    nlev=n_elements(tt)/nlon/nlat
    tt=reform(tt, nlon*nlat, nlev, /overwrite)
    pp=reform(pp, nlon*nlat, nlev, /overwrite)
    if n_elements(s1) ne 1 then ecmwf_iso_settings
    ;t0=tt[*, t0lev]
    tlapse=fltarr(nlon*nlat)
    min_corr=1.
    t0=fltarr(nlon*nlat)
    for i=0L, nlon*nlat-1 do begin
      tprof=reform(tt[i, s1:s2])
      pprof=reform(pp[i, s1:s2])
      tlapse[i]=regress(alog(tprof), alog(pprof))
      tlapse[i]=tlapse[i]/9.8
      t0[i]=total(tprof)/n_elements(tprof)
      corr=correlate(alog(tprof), alog(pprof))
      if corr lt min_corr then min_corr=corr
    endfor
    print, "Minimum correlation between p and t:", min_corr

    ;generate the retrieval matrix:
    if file_test(base_name+"chan18_bt"+datestr, /regular, /read) ne 1 then ecmwf_field_amsu_bt, date
    chan18=fltarr(nlon*nlat)
    openr, lun, base_name+"chan18_bt"+datestr, /get_lun
    readu, lun, chan18
    free_lun, lun
    chan19=fltarr(nlon*nlat)
    openr, lun, base_name+"chan19_bt"+datestr, /get_lun
    readu, lun, chan19
    free_lun, lun
    chan20=fltarr(nlon*nlat)
    openr, lun, base_name+"chan20_bt"+datestr, /get_lun
    readu, lun, chan20
    free_lun, lun
    chan21=fltarr(nlon*nlat)
    openr, lun, base_name+"chan21_bt"+datestr, /get_lun
    readu, lun, chan21
    free_lun, lun
    chan22=fltarr(nlon*nlat)
    openr, lun, base_name+"chan22_bt"+datestr, /get_lun
    readu, lun, chan22
    free_lun, lun

    retmat=transpose([[t0], [tlapse], [chan18], [chan19], [chan20], [chan21], [chan22]])
    retmat=retmat[variables, *]

    ;store the retrieval matrix:
    openw, lun, retbase+".vec", /get_lun
    writeu, lun, n_elements(variables)
    writeu, lun, retmat
    free_lun, lun
  endif

  ;get the training data:
  if file_test(trainbase+".vec", /read, /regular) ne 1 then ecmwf_training_matrix, date
;  openr, lun, "train_mat"+datestr+".dat", /get_lun
;  openr, lun, filename, /get_lun
;  nvar=0L
;  readu, lun, nvar
;  ntrain=((fstat(lun)).size-4)/nvar/4
;  train_mat=fltarr(nvar, ntrain)
;  readu, lun, train_mat
;  free_lun, lun

  ;train_mat=train_mat[[0, 2, 3, 4, 5, 6], *]

  openr, lun, base_name+"vmr_train"+sstr+datestr+".dat", /get_lun
  ntrain=(fstat(lun)).size/4
  vmr=fltarr(ntrain)
  readu, lun, vmr
  free_lun, lun

  classes=vmr gt vmr_thresh

  ;write out the classes:
  tstr=string(vmr_thresh)
  tstr=strtrim(tstr, 2)
  openw, lun, trainbase+".cls"
  writeu, lun, long(classes)
  free_lun, lun

  stop

  outbase=retbase

  if file_test(outbase+".cls", /read, /regular) ne 1 then begin
    print, "Starting the classifier..."
    spawn, "time classify_a-O2 "+trainbase+" "+retbase+".vec "+outbase+string(knn)+string(wopt)
  endif

  ;classify, train_mat, retmat, classes, 2, field, confidence, width=width
  ;classify_caller, train_mat[varind, *], retmat[varind, *], classes, knn, wopt, field, confidence

  field=lonarr(nlon, nlat)
  confidence=fltarr(nlon, nlat)

  openr, lun, outbase+".cls", /get_lun
  readu, lun, field
  free_lun, lun
  openr, lun, outbase+".con", /get_lun
  readu, lun, confidence
  free_lun, lun

end

pro ecmwf_iso_plot_comp, date, field, sigma, vmr_thresh, confidence, $
			psfile=psfile, sat=sat, tiffile=tiffile
  get_ecmwf_field, date, "SQ", sh, lon, lat, $
		path="/home/millsp/data/ecmwf_misc"
  vmr=sh[*, *, sigma];*1.60627

  ;write this to a file for later use:
  time=convert_date(date)
  datestr=string(time.year, time.month, time.day, time.hour, $
  		format="(i4.4, i2.2, i2.2, i2.2)")
  sstr=string(sigma, format="(i2.2)")+"s"
  openw, lun, "vmr"+sstr+datestr, /get_lun
  writeu, lun, vmr
  free_lun, lun

  nlon=n_elements(lon)
  nlat=n_elements(lat)
  lon=rebin(lon, nlon, nlat)
  lat=rebin(transpose(lat), nlon, nlat)

  truth=vmr gt vmr_thresh
  r=correlate(truth, field)
  print, "r=", r
  acc=total(truth eq field)/nlon/nlat

  title="ECMWF"+string(vmr_thresh)+" vmr isoline for "+time_string(convert_date(date))+" sigma level"+string(sigma)
;  subtitle="retrieved (dots) vs. actual (thick solid line); r="+string(r)
  subtitle="r="+string(r)+";  acc.:"+string(acc*100)+"%"

;  map_set, 0, 0, 0, /cylind, title=title
;  map_continents
;  map_grid, londel=30, latdel=30, glinestyle=1

  theta=findgen(40)*2*!pi/40
  symx=cos(theta)
  symy=sin(theta)
  usersym, symx, symy

  ;create the levels:
  nlev=20
  conlev=(findgen(nlev))/(nlev-1)

  ;create the color index:
  ncol=nlev
  ;green=[0, fltarr(ncol), 250, 200, 200]
  ;red=[findgen(ncol+1)*200/(ncol), 250, 0, 221]
  ;blue=[0, 200-findgen(ncol)*200/(ncol-1), 250, 0, 50]
  
  red=round([0, findgen(nlev)/(nlev-1)*155+50, 255, 200, 200])
  blue=round([0, 255-findgen(nlev)/(nlev-1)*255, 255, 0, 221])
  green=round([0, 150-findgen(nlev)/(nlev-1)*100, 255, 0, 50])
  
  colind=rgbtoindex(red, green, blue)

  ;ret=rgbtoindex(221, 200, 50)
  ret=rgbtoindex(0, 0, 0)

  if n_elements(psfile) eq 1 then begin
    set_plot, "ps"
    device, filename=psfile
    device, /color, bits_per_pixel=8
    device, /landscape, xsize=9, ysize=6.5, xoffset=1, yoffset=10, /inches
    tvlct, fix(red), fix(green), fix(blue), 1
    colind=indgen(nlev+3)+1
  endif
  
  white=colind[ncol+1]
  black=colind[0]
  ccol=colind[1:nlev]

  if n_elements(tiffile) eq 1 then window, 0, xsize=900, ysize=600

  plot, lon, lat, title=title, subtitle=subtitle, xtitle="lon (deg.)", ytitle="lat (deg.)", $
  		/nodata, xstyle=1, ystyle=1, background=white, color=black, charsize=1.3, $
		charthick=2

  plot_continents, color=black

  if n_elements(confidence) eq nlon*nlat then begin

    contour, confidence, lon, lat, levels=conlev, c_color=ccol, /cell_fill, /overplot, $
	    	background=white, color=black

    ;plot each point with a different color based on the confidence rating:
    ;for i=0, nlon*nlat-1 do begin
    ;  if field[i] eq 1 then begin
    ;    if finite(confidence) then begin
;	  color=colind[confidence[i]*(ncol-1)]
;	endif else begin
;	  color=0
;	endelse
    ;    plots, lon[i], lat[i], color=color, psym=2, symsize=0.15, /data
    ;  endif
    ;endfor
    plot_continents, color=black
  endif

  if keyword_set(sat) then begin
    get_ecmwf_field, date, "TT", tt, lon, lat
    get_ecmwf_field, date, "PP", pp, lon, lat
    ;calculate dewpoint:
    t=tt[*, *, sigma]
    p=pp[*, *, sigma]
    e=p
    indi=where(t lt 273.15, cnti)
    indw=where(t ge 273.15, cntw)
    if cnti ne 0 then e[indi]=eq_vp_ice(t[indi])
    if cntw ne 0 then e[indw]=eq_vp(t[indw])
    wvp=vmr*pp[*, *, sigma]*100
    satind=reform(wvp ge e, nlon, nlat)
    contour, satind, lon, lat, levels=1, c_color=colind[nlev+2], /fill, $
	    	background=white, color=green
  endif

  levels=[0.6, 0.8, 1, 1.2, 1.4]*vmr_thresh
  contour, vmr, lon, lat, levels=levels, /overplot, c_thick=[1, 1, 5, 1, 1], $
	  	c_color=ret

  ind=where(field eq 1, cnt)
  if cnt ne 0 then begin
    x=lon[ind]
    y=lat[ind]
    oplot, x, y, psym=8, symsize=0.15, color=ret
  endif

  if n_elements(psfile) eq 1 then begin
    device, /close
    set_plot, "x"
  endif

  if n_elements(tiffile) eq 1 then begin
    image=tvrd(0)
    ind=where(image eq 250)
    indr=where(image eq 221)
    image=fix(round(image*(nlev-1.)/200))
    image[ind]=nlev+1
    image[indr]=nlev+3

    write_tiff, tiffile, image, $
              red=fix(red), blue=fix(blue), green=fix(green)
  endif

  stop

end

pro ecmwf_iso_plot_range, date, field, sigma, vmr_thresh, confidence, min_con, $
			psfile=psfile, sat=sat, tiffile=tiffile
  get_ecmwf_field, date, "SQ", sh, lon, lat, $
		path="/home/millsp/data/ecmwf_misc"
  vmr=sh[*, *, sigma];*1.60627

  ;write this to a file for later use:
  time=convert_date(date)
  datestr=string(time.year, time.month, time.day, time.hour, $
  		format="(i4.4, i2.2, i2.2, i2.2)")
  sstr=string(sigma, format="(i2.2)")+"s"
  openw, lun, "vmr"+sstr+datestr, /get_lun
  writeu, lun, vmr
  free_lun, lun

  nlon=n_elements(lon)
  nlat=n_elements(lat)
  lon=rebin(lon, nlon, nlat)
  lat=rebin(transpose(lat), nlon, nlat)

  truth=vmr gt vmr_thresh
  r=correlate(truth, field)
  print, "r=", r

  title="ECMWF"+string(vmr_thresh)+" vmr isoline for "+time_string(convert_date(date))+" sigma level"+string(sigma)
  ;subtitle="retrieved (dots) vs. actual (thick solid line); r="+string(r)

;  map_set, 0, 0, 0, /cylind, title=title
;  map_continents
;  map_grid, londel=30, latdel=30, glinestyle=1

  white=rgbtoindex(250, 250, 250)
  black=rgbtoindex(0, 0, 0)
  red=rgbtoindex(200, 0, 0)
  grey=rgbtoindex(100, 100, 100)

  if n_elements(psfile) eq 1 then begin
    set_plot, "ps"
    device, filename=psfile
    device, /color, bits_per_pixel=8
    device, /landscape, xsize=9, ysize=6.5, xoffset=1, yoffset=10, /inches
    tvlct, fix(red), fix(green), fix(blue), 1
    colind=indgen(nlev+1)+1
    ccol=colind[1:nlev]
  endif

  if n_elements(tiffile) eq 1 then window, 0, xsize=900, ysize=600

  inner_bound=field and confidence gt min_con
  outer_bound=field or confidence lt min_con
  range=outer_bound and inner_bound ne 1

;  contour, range, lon, lat, levels=0.5, c_color=grey, /fill, /overplot, $
;  		background=white

  contour, confidence, lon, lat, levels=[0, min_con], c_color=[grey, white], $
  		title=title, subtitle=subtitle, xtitle="lon (deg.)", ytitle="lat (deg.)", $
  		xstyle=1, ystyle=1, background=white, color=black, charsize=1.5, $
		charthick=2, /cell_fill

  plot_continents, color=black

  contour, confidence, lon, lat, levels=[0, min_con], c_linestyle=[2, 2], $
	  	/overplot, color=black

;  contour, inner_bound, lon, lat, levels=0.5, /overplot, $
;	  		color=black, c_linestyle=2, c_thick=1
;  contour, outer_bound, lon, lat, levels=0.5, /overplot, $
;	  		color=black, c_linestyle=2, c_thick=1
  contour, field, lon, lat, levels=0.5, color=black, c_thick=2, /overplot

  levels=[0.6, 0.8, 1, 1.2, 1.4]*vmr_thresh
  contour, vmr, lon, lat, levels=vmr_thresh, /overplot, $
	  	color=red, c_thick=2

  if n_elements(psfile) eq 1 then begin
    device, /close
    set_plot, "x"
  endif

  if n_elements(tiffile) eq 1 then begin
    image=tvrd(0)
    indr=where(image eq 200)
    indb=where(image eq 0)
    indw=where(image eq 250)
    indg=where(image eq 100)
    image[indb]=0
    image[indw]=1
    image[indr]=2
    image[indg]=3

    write_tiff, tiffile, image, $
              red=[0, 250, 200, 100], blue=[0, 250, 0, 100], green=[0, 250, 0, 100]
  endif

  stop

end

pro ecmwf_iso_con_stats, date, sigma, vmr_thresh, field, confidence

  white=rgbtoindex(250, 250, 250)
  black=rgbtoindex(0, 0, 0)

  get_ecmwf_field, date, "SQ", sh, lon, lat
  vmr=sh[*, *, sigma]*1.60627

  truth=vmr gt vmr_thresh

  ;calculate histograms for number of correct classifications:
  hmin=0
  hmax=1
  nbins=21
  bins=findgen(nbins-1)/20+0.025

  ind=where(finite(confidence))
  con1=confidence[ind]
  truth1=truth[ind]
  field1=field[ind]

  tind=where(truth1 eq field1, nt)
  find=where(truth1 ne field1, nf)

  thist=histogram(con1(tind), min=hmin, max=hmax, nbins=nbins)
  fhist=histogram(con1(find), min=hmin, max=hmax, nbins=nbins)

;  window, 0

;  plot, bins, thist
;  oplot, bins, fhist

  ;window, 1, xsize=800, ysize=550

  set_plot, "ps"
  device, filename="acc_vs_con3.ps"
  device, /landscape, xsize=9, ysize=6.5, xoffset=1, yoffset=10, /inches

  plot, bins, 1.0*thist/(thist+fhist), psym=2, xtitle="confidence rating", ytitle="accuracy", $
  		title="Classification accuracy vs. confidence rating", /ynozero, charsize=1.5, $
		charthick=1.5, thick=1.5
  oplot, [0, 1], [0.5, 1], thick=1.5
  axis, 1, 0, 0, /yaxis, yrange=[0, 1], yminor=10, charsize=1.5

;  image=tvrd(0)
;  indb=where(image eq 0)
;  indw=where(image eq 255)
;  image[indb]=0
;  image[indw]=1
;  write_tiff, "acc_v_conr.tif", image, red=[50, 250], blue=[50, 250], green=[50, 250]

  tsi=sort(con1[tind])
  fsi=sort(con1[find])

  tcon=con1[tind[tsi]]
  fcon=con1[find[fsi]]

;  window, 2

;  plot, tcon, nt-findgen(nt)
;  oplot, fcon, nf-findgen(nf)

  ;window, 3, xsize=800, ysize=550

  congrd=findgen(101)/100
  pt=interpol(nt-findgen(nt), tcon, congrd)
  pf=interpol(nf-findgen(nf), fcon, congrd)
  pt[100]=1
  pf[100]=1

  ;approximate formulae:
  conint1=1-1.0*pf/(nt-pt+pf)
  conint2=1-1.0*pf/(nt-pt-nf+2*pf)
  conint3=1-1.0*pf/(nf)
  conint10=congrd*conint2+(1-congrd)*conint3

  ;compare with a tedious (but more accurate) line integral:
  ecmwf_iso_line_int, vmr, vmr_thresh, confidence, path_cons, conint4

  set_plot, "ps"
  device, /close_file
  device, filename="con_int_func3.ps"
  device, /landscape, xsize=9, ysize=6.5, xoffset=1, yoffset=10, /inches

  plot, path_cons, conint4, xtitle="confidence rating", $
	  	ytitle="confidence interval", yrange=[0, 1], $
		title="Isoline vs. classification confidence", charsize=1.5, thick=1.5
;  oplot, congrd, conint2
;  oplot, congrd, conint3
;  oplot, path_cons, conint4
;  oplot, congrd, conint10, linestyle=2

  c90=interpol(path_cons, conint4, 0.9)
  c50=interpol(path_cons, conint4, 0.5)
  cstd=interpol(path_cons, conint4, 0.68269)

  oplot, [0, c90], [0.9, 0.9]
  oplot, [c90, c90], [0, 0.9]
  oplot, [0, c50], [0.5, 0.5]
  oplot, [c50, c50], [0, 0.5]
  oplot, [0, cstd], [0.68269, 0.68269]
  oplot, [cstd, cstd], [0, 0.68260]

  xyouts, 0.05, 0.91, '90 %', /data, charsize=1.3
  xyouts, 0.05, 0.51, '50 %', /data, charsize=1.3
  xyouts, 0.05, 0.69269, '1 std. dev.', /data, charsize=1.3

;  image=tvrd(3)
;  indb=where(image eq 0)
;  indw=where(image eq 255)
;  image[indb]=0
;  image[indw]=1
;  write_tiff, "iso_v_class_con.tif", image, red=[50, 250], blue=[50, 250], green=[50, 250]

;device, /close
;set_plot, "x"

  device, /close_file
  device, filename="isoline_con3.ps"
  device, /color, bits_per_pixel=8

  tvlct, [250, 150], [0, 150], [0, 150], 1
  red=1
  grey=2

  ;window, 4, xsize=900, ysize=600
;  levels=[0.6, 0.8, 1, 1.2, 1.4]*vmr_thresh
  levels=vmr_thresh
  map_set, 0, 180, 0, /cylind, title="Isoline confidence interval"
  map_continents, color=grey, /fill
;  contour, vmr, lon, lat, levels=levels, /nodata, color=white, $
;	  	title="Isoline confidence interval", charsize=1.5, $
;	  	background=rgbtoindex(201, 201, 201), $
;	  	xstyle=1, ystyle=1, xtitle="lon [deg.]", ytitle="lat [deg.]"
;  plot_continents, color=grey, /fill
  contour, vmr, lon, lat, levels=levels, color=red, $
  		c_thick=5, /overplot
  contour, confidence, lon, lat, levels=[c50, cstd, c90], c_linestyle=[0, 2, 0], $
  		/overplot;, color=black
  contour, field, lon, lat, levels=0.5, c_thick=5, /overplot, c_color=black
;  contour, vmr, lon, lat, levels=levels, color=rgbtoindex(200, 0, 0), $
;  		c_thick=[1, 1, 2, 1, 1], /overplot
;  image=tvrd(4)
;  indw=where(image eq 250)
;  indb=where(image eq 0)
;  indg=where(image eq 201)
;  indr=where(image eq 251)
;  image[indw]=0
;  image[indb]=1
;  image[indg]=2
;  image[indr]=3

;  write_tiff, "iso_con_map.tif", image, red=[250, 0, 201, 250], $
;	  	blue=[250, 0, 201, 50], green=[250, 0, 201, 50]

  device, /close
  set_plot, "x"

  stop

end

;do a tedious line-integral to get the confidence interval:
pro ecmwf_iso_line_int, vmr, vmr_thresh, confidence, path_cons, conint

  contour, vmr, levels=vmr_thresh, /path_data_coords, path_xy=path_xy, $
	  	path_info=path_info

  npts=n_elements(path_xy)/2

  path_con=interpolate(confidence, path_xy[0, *], path_xy[1, *])

  ninfo=n_elements(path_info)

  ds=fltarr(npts)
  for i=0L, ninfo-1 do begin
    j1=path_info[i].offset
    n=path_info[i].n-1
    j2=j1+n-1
    ds1=point_spacing(path_xy[0, j1:j2], path_xy[1, j1:j2], /nowrap)
    if n gt 2 then begin
      ds[j1]=[ds1[0]/2, (ds1[0:n-3]+ds1[1:n-2])/2, ds1[n-2]/2]
    endif else begin
      ds[j1]=[ds1[0]/2, ds1[0]/2]
    endelse
  endfor

  s=total(ds)
  sind=sort(path_con)
  dss=ds[sind]
  
  conint=fltarr(npts)
  conint[0]=dss[0]
  for i=1L, npts-1 do conint[i]=conint[i-1]+dss[i]
  conint=conint/s

  path_cons=path_con[sind]

end

