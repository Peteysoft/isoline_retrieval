@~/climate/block2/file_ut

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
   
;program to animate a years worth of data in five-day averages           
           
pro plt_N, grid, longitude, latitude, levels=contour_levels, $           
           fill=fill_values, n_levels=n_levels           
  ;produces a map of the Northern hemisphere with a contour plot overlaid           
  ;"stripped down" version of the "plot_north" procedure           
             
  ;print, grid           
           
  if n_elements(contour_levels) eq 0 then begin           
    max_val=max(grid)           
    min_val=min(grid)           
    if n_elements(n_levels) eq 0 then n_levels=11           
    contour_levels=min_val+findgen(n_levels)*((max_val-min_val)/(n_levels-1))           
  endif else n_levels=n_elements(contour_levels)           
               
  ;plot map, contour, and legend           
           
  map_set, 90, 0, 0, /stereo, limit= [0,-180,90,180]           
             
;  print, longitude           
;  print, latitude           
             
  contour, grid, longitude, latitude, /overplot, $           
        /fill, c_colors=fill_values, levels=contour_levels           
           
  map_continents;, /hires           
           
end           
;****************************************************************************           
           
pro animate, file, start, finish, inter=inter, height=height, width=width, $   
           ratio=ratio, levels=levels
             
  if n_elements(ratio) eq 0 then ratio=1   
           
  ;next read the data from a file:           
  data=open_assoc(file, handle, n_times=n_times)           
  lon=dim_data(file, 0)           
  lat=dim_data(file, 1)
  east=where(lon lt 180)
  west=where(lon ge 180)
  if west(0) ne -1 then begin
    lon_fixed=[lon(west)-360, lon(east)] 
  endif else lon_fixed=lon
  lon_range=[lon_fixed(0), lon_fixed(n_elements(lon)-1)]   
  lat_range=[lat(0), lat(n_elements(lat)-1)]   
  print, lon_range   
  print, lat_range   
  if n_elements(width) eq 0 then width=300   
  if n_elements(height) eq 0 then begin   
    height=width*(lat_range(1)-lat_range(0))/(lon_range(1)-lon_range(0))   
  endif   
            
  if n_elements(finish) eq 0 then begin        
    if n_elements(start) ne 0 then finish=start else finish=n_times-1        
  endif        
  if n_elements(start) eq 0 then start=0           
           
  ;load colour table and set fill values:           
  n_levels=11           
  device, decomposed=0           
  loadct, 33, bottom=1, ncolors=n_levels           
  fill_values=indgen(n_levels)+1           
                           
  ;now to plot all the images           
  n_maps=(finish-start)/ratio     
  window, 1, xsize=width, ysize=height, xpos=200, ypos=0   
  if keyword_set(inter) then begin             
    xinteranimate, set=[width,height,n_maps+1]  
  endif           
  for t=0, n_maps do begin
    if west(0) ne -1 then begin
      map=[data(west, *, start+t*ratio), data(east, *, start+t*ratio)]
    endif else map=data(start+t*ratio)
    if n_elements(levels) eq 0 then begin           
      max_val=max(map)           
      min_val=min(map)           
      levels=min_val+findgen(n_levels)*((max_val-min_val)/(n_levels-1))           
    endif       
    contour, map, lon_fixed, lat, levels=levels, /fill, c_colors=fill_values, $     
           xmargin=[0,0], ymargin=[0,0], xstyle=4, ystyle=4      
    plot_continents, lon_range=lon_range, lat_range=lat_range   
;    plt_N, map, lon, lat, levels=levels, fill=fill_values, n_levels=n_levels           
    if keyword_set(inter) then xinteranimate, frame=t, window=1           
  endfor           
             
  ;finally display the animated images           
  if keyword_set(inter) then xinteranimate           
          
  free_lun, handle          
             
end           


           
