;@~/climate/block2/file_ut

;the idl version doesn't work, so I've taken it upon myself to write   
;my own:   
pro plot_continents, lon_range=lon_range, lat_range=lat_range, color=color   
  if n_elements(lon_range) eq 0 then lon_range=!x.crange
  if n_elements(lat_range) eq 0 then lat_range=!y.crange
  n_segments=283   
  mapfilename="/opt/rsi/idl/resource/maps/supmap.dat"   
  if lon_range[0] ge 0. and lon_range[1] gt 180 then lon_corr=1 else lon_corr=0
   
  openr, mapfile, mapfilename, /xdr, /get_lun   
   
  for segment=0, n_segments do begin   
    npts=0L   
    maxlat=0.0   
    minlat=0.0   
    maxlon=0.0   
    minlon=0.0   
    readu, mapfile, npts
    readu, mapfile, maxlat
    readu, mapfile, minlat
    readu, mapfile, maxlon
    readu, mapfile, minlon 
    if lon_corr then begin
      if maxlon lt 0 then maxlon=maxlon+360
      if minlon lt 0 then minlon=minlon+360
    endif
    if minlon gt lon_range(1) or maxlon lt lon_range(0) or $   
       minlat gt lat_range(1) or maxlat lt lat_range(0) then begin   
      point_lun, -mapfile, current   
      point_lun, mapfile, current+npts*4   
    endif else begin   
      outline=replicate({lat:0.0, lon:0.0}, npts/2)
      readu, mapfile, outline
      lon=outline.lon
      if lon_corr then begin
        ind1=where(lon lt 0, cnt1, complement=ind2, ncomplement=cnt2)
	if cnt1 ne 0 then oplot, 360+lon[ind1], outline.lat, color=color
	if cnt2 ne 0 then oplot, lon[ind2], outline.lat, color=color
      endif else begin
        oplot, lon, outline.lat, color=color
      endelse
    endelse   
  endfor   
  free_lun, mapfile   
end   
   
