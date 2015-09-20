@wv_formulas
@time_lib
;@animate
@boundary_geometry

pro ecmwf_iso_plot_comp, date, field, sigma, vmr_thresh, confidence, $
			psfile=psfile, sat=sat, tiffile=tiffile
  get_ecmwf_field, date, "SQ", sh, lon, lat
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

  title="ECMWF"+string(vmr_thresh, format="(g8.2)")+" mmr isoline for "+time_string(convert_date(date))+" sigma level"+string(sigma, format="(i3)")
;  subtitle="retrieved (dots) vs. actual (thick solid line); r="+string(r)
;  subtitle="r="+string(r)+";  acc.:"+string(acc*100)+"%"
  subtitle="accuracy = "+string(acc*100, format="(f6.2)")+"%"

;  map_set, 0, 0, 0, /cylind, title=title
;  map_continents
;  map_grid, londel=30, latdel=30, glinestyle=1

  theta=findgen(40)*2*!pi/40
  symx=cos(theta)
  symy=sin(theta)
  usersym, symx, symy, /fill

  ;create the levels:
  nlev=21
  conlev=(findgen(nlev))/(nlev-1)

  ;create the color index:
  ncol=nlev
  ;green=[0, fltarr(ncol), 250, 200, 200]
  ;red=[findgen(ncol+1)*200/(ncol), 250, 0, 221]
  ;blue=[0, 200-findgen(ncol)*200/(ncol-1), 250, 0, 50]
  
  ;red=round([0, findgen(nlev)/(nlev-1)*155+50, 255, 200, 200])
  ;blue=round([0, 255-findgen(nlev)/(nlev-1)*255, 255, 0, 221])
  ;green=round([0, 150-findgen(nlev)/(nlev-1)*100, 255, 0, 50])
 
  red=findgen(nlev)*255/(nlev+1)
  green=fltarr(nlev)
  blue=255-red
  
  colind=rgbtoindex(red, green, blue)

  ;ret=rgbtoindex(221, 200, 50)
  ret=rgbtoindex(0, 0, 0)

  if n_elements(psfile) eq 1 then begin
    set_plot, "ps"
    device, filename=psfile
    device, /color, bits_per_pixel=8
    device, /landscape, xsize=9, ysize=6.5, xoffset=1, yoffset=10, /inches
    device, /helvetica
    tvlct, fix(red), fix(green), fix(blue), 1
    colind=indgen(nlev+1)
  endif
  
  ;white=colind[ncol+1]
  ;black=colind[0]
  ccol=colind[1:nlev]

  if n_elements(tiffile) eq 1 then window, 0, xsize=900, ysize=600

  plot, lon, lat, title=title, subtitle=subtitle, xtitle="lon (deg.)", ytitle="lat (deg.)", $
  		/nodata, xstyle=1, ystyle=1, background=white, color=black, charsize=1.3, $
		charthick=3

  ;plot_continents, color=black

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
  endif

  plot_continents, color=255

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
    contour, satind, lon, lat, levels=1, c_color=colind[nlev+1], /fill, $
	    	background=white, color=green
  endif

  levels=[0.6, 0.8, 1, 1.2, 1.4]*vmr_thresh
  contour, vmr, lon, lat, levels=levels, /overplot, c_thick=[1, 1, 5, 1, 1];, $
;	  	c_color=255

  ind=where(field eq 1, cnt)
  if cnt ne 0 then begin
    x=lon[ind]
    y=lat[ind]
    oplot, x, y, psym=8, symsize=0.2;, color=255
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
  get_ecmwf_field, date, "SQ", sh, lon, lat
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

  title="ECMWF"+string(vmr_thresh, format="(g8.2)")+" mmr isoline for "+time_string(convert_date(date))+" sigma level"+string(sigma, format="(i3)")
  ;subtitle="retrieved (dots) vs. actual (thick solid line); r="+string(r)

;  map_set, 0, 0, 0, /cylind, title=title
;  map_continents
;  map_grid, londel=30, latdel=30, glinestyle=1

;  white=rgbtoindex(250, 250, 250)
  white=255
  ;black=rgbtoindex(0, 0, 0)
  red=rgbtoindex(200, 0, 0)
  grey=rgbtoindex(100, 100, 100)

  if n_elements(psfile) eq 1 then begin
    set_plot, "ps"
    device, filename=psfile
    ;device, /color, bits_per_pixel=8
    device, /landscape, xsize=9, ysize=6.5, xoffset=1, yoffset=10, /inches
    device, /helvetica
    ;tvlct, fix(red), fix(green), fix(blue), 1
    ;colind=indgen(nlev+1)+1
    ;ccol=colind[1:nlev]
    grey=150
  endif

  if n_elements(tiffile) eq 1 then window, 0, xsize=900, ysize=600

  inner_bound=field and confidence gt min_con
  outer_bound=field or confidence lt min_con
  range=outer_bound and inner_bound ne 1

;  contour, range, lon, lat, levels=0.5, c_color=grey, /fill, /overplot, $
;  		background=white

  contour, confidence, lon, lat, levels=[0, min_con], c_color=[grey, white], $
  		title=title, subtitle=subtitle, xtitle="lon (deg.)", ytitle="lat (deg.)", $
  		xstyle=1, ystyle=1, charsize=1.3, $
		charthick=3, /cell_fill;,
;		background=white, color=black

  plot_continents;, color=black

;  contour, confidence, lon, lat, levels=[0, min_con], c_linestyle=[2, 2], $
;	  	/overplot, color=black

;  contour, inner_bound, lon, lat, levels=0.5, /overplot, $
;	  		color=black, c_linestyle=2, c_thick=1
;  contour, outer_bound, lon, lat, levels=0.5, /overplot, $
;	  		color=black, c_linestyle=2, c_thick=1
;  contour, field, lon, lat, levels=0.5, color=black, c_thick=2, /overplot
  contour, field, lon, lat, levels=0.5, c_linestyle=1, c_thick=6, /overplot

;  levels=[0.6, 0.8, 1, 1.2, 1.4]*vmr_thresh
  contour, vmr, lon, lat, levels=vmr_thresh, /overplot, $
	  	c_thick=4

  oplot, [20, 40], [-55, -55], thick=4
  oplot, [20, 40], [-65, -65], linestyle=1, thick=6
  polyfill, [20, 40, 40, 20], [-72, -72, -77, -77], color=grey

  xyouts, 50, -55, "actual", charsize=1.3, /data, charthick=2
  xyouts, 50, -65, "retrieved", charsize=1.3, /data, charthick=2
  xyouts, 50, -75, "90% tolerance", charsize=1.3, /data, charthick=2

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


