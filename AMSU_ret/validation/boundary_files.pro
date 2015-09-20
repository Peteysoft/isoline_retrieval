;a set of utilities for working with "boundary" files:


;scans a "boundary file": returns an array of time values, an array of sizes
;and an array of locations:
function scan_boundary_file, filename, nall, loc=loc

;+**************************************************************************
;		SCAN_BOUNDARY_FILE
;***************************************************************************
;
; usage:	time=scan_boundary_file(filename [, elements][, loc=loc])
;
; purpose:	Scans a binary file containing a series of evected contours
;		or "boundaries."  Returns an array of time values (see time_lib.pro)
;		for each record.  Also returns the number of elements in each
;		record and the location of the start of each record.  Returns
;		an error code if the routine fails.
;
; parameters:	filename	The name of the file to scan.
;
;		elements	An array containing the number of longitude and latitude 
;					values in each record
;
; keywords:	loc			An array containing the location of the start of each 
;					record in the file.  Each record has seven headers.  
;					The first six headers represent the time for which the 
;					boundary is valid and are organized as follows: year, month, 
;					day, hour, minute and seconds.  All of them are represented 
;					as two byte integers except for the seconds field which
;					is a four byte floating point value.  The last header
;					is the number of longitude and latitude values and is
;					stored as a four-byte integer.  All the 
;					longitude values and all the latitude values are stored 
;					contiguously, NOT as lon-lat pairs.  The files themselves
;					have two headers, both longword integers: a "magic"
;					number which is set at -55 and the number of records
;					in the file.  The potential temperature level is not stored.
;					(note to self... )
;
; dependencies:	time_lib.pro
;
; history:	First formally documented 2002-01-29 Peter Mills (peter.mills@nrl.navy.mil)
;
;-***************************************************************************

  ;open the file:
  openr, lun, filename, /get_lun, err=err
  if err ne 0 then begin
    message, 'File "'+filename+'" not found', /continue
    return, err
  endif
  
  on_ioerror, fail

  ;get the "magic" number:
  magic=0L & nt=0L
  readu, lun, magic
  if magic ne -55L then begin
    message, "Likely file error, magic number incorrect--"+ $
    		string(magic), /continue
  endif
  
  ;get the number of records and allocate the arrays:
  readu, lun, nt
  alltimes=replicate({time_str}, nt)
  nall=lonarr(nt)
  loc=lonarr(nt)
  
  ;scan through each record in the file, storing the appropriate values:
  time={time_str}
  n=0L
  newloc=8
  for i=0L, nt-1 do begin
    if eof(lun) then begin
      message, string(i+1)+" records found; does not agree with file header:"+ $
        		string(nt), /continue
      alltimes=alltimes[0:i]
      nall=nall[0:i]
      loc=loc[0:i]
      break
    endif
    loc[i]=newloc
    readu, lun, time, n
    alltimes[i]=time
    nall[i]=n
    point_lun, -lun, pos
    newloc=pos+8*n
    point_lun, lun, newloc
  endfor
  
  free_lun, lun
  
  return, alltimes
  
  fail:
    free_lun, lun
    message, "I/O error", /continue
    return, !error_state.code
  
end

;gets the nth boundary from a "boundary evolution" file:
;if n=-1 then gets the boundary defined at or later than time
pro get_boundary, file, n0, time, lon, lat, err=err

;+**************************************************************************
;		GET_BOUNDARY
;***************************************************************************
;
; usage:	get_boundary, filename, index, time, lon, lat, err=err
;
; purpose:	Scans a binary file containing a series of evected contours
;		or "boundaries" for a specific record.  Returns the two arrays
;		of longitude and latitude values as well as the time for which
;		the boundary is valid.
;
; parameters:	filename	The name of the file to read.
;
;		index		The zero-based index of the record to retrieve from
;					the file.  If this is a -1, then searches the file
;					for a time instead and stores the index of the record
;					actually retrieved in this parameter.
;
;		time		The time for which the boundary is valid.  If index
;					is -1 then this is taken as an input parameter.  The
;					routine then searches the file for the first record
;					valid on or later than that time.  Time values are
;					represented as structures with six fields; see the
;					file "time_lib.pro"
;
;		lon			A returned array of longitude values.
;
;		lat			A returned array of latitude values.
;
; dependencies:	time_lib.pro
;
; history:	First formally documented 2002-01-29 Peter Mills (peter.mills@nrl.navy.mil)
;
;-***************************************************************************

  openr, lun, file, /get_lun, err=err
  if err ne 0 then begin
    message, 'File "'+file+'" not found', /continue
    return
  endif
  
  on_ioerror, fail

  magic=0L & nt=0L
  readu, lun, magic, nt
  if magic ne -55L then begin
    message, "Likely file error, magic number incorrect--"+ $
    		string(magic), /continue
  endif

  n=0L
  
  if n0 eq -1 then begin
    time1={time_str}
    readu, lun, time1, n
    while time_compare(time1, time) lt 1 do begin
      point_lun, -lun, pos
      point_lun, lun, pos+8*n
      if eof(lun) then begin
        message, "Time value out of bounds--"+ $
        		time_string(time), /continue
        free_lun, lun
        err=-1
        return
      endif
      readu, lun, time1, n
      n0=n0+1
    endwhile
    lon=fltarr(n)
    lat=fltarr(n)
    readu, lun, lon, lat
    n0=n0+1
    time=time1
  endif else begin    
  
    time={time_str}
    for i=0, n0-1 do begin
      readu, lun, time, n
      point_lun, -lun, pos
      point_lun, lun, pos+8*n
      if eof(lun) then begin
        message, "Index value out of bounds--"+ $
        		string(n0), /continue
        free_lun, lun
        err=-1
        return
      endif
    endfor
    readu, lun, time, n
    lon=fltarr(n)
    lat=fltarr(n)
    readu, lun, lon, lat
  endelse
  
  free_lun, lun
  err=0
  return
  
  fail:
    free_lun, lun
    message, "I/O error", /continue
    err=!error_state.code
    return

end


;if the number of saved boundaries at the beginning of the file is wrong,
;fixes it:
pro fix_boundary_file, file
  openu, lun, file, /get_lun
  
  magic=0L & nt=0L
  readu, lun, magic, nt

  on_ioerror, fixnow

  time={time_str}
  n=0L
  i=0L
  while not eof(lun) do begin
    readu, lun, time, n
    lon=fltarr(n)
    lat=fltarr(n)
    readu, lun, lon, lat
    i=i+1
;    stop
  endwhile
  
;  free_lun, lun
  
;  openu, lun, file, /get_lun

  fixnow:

  print, i, " records found in file: ", file
  
  point_lun, lun, 4
  writeu, lun, i
  
  free_lun, lun

end

;adds in "magic" number and number of time steps to a file containing
;boundary evolution
pro convert_boundary_file, oldfile, newfile
  magic=-55L

  openr, inlun, oldfile
  openw, outlun, newfile
  writeu, outlun, magic
  ;write a "place-holder" for the number of time-steps:
  writeu, outlun, 0L
  
  time={time_str}
  n=0
  readu, inlun, time, n
  lon=fltarr(n)
  lat=fltarr(n)
  readu, inlun, lon, lat
  writeu, outlun, time, long(n), lon, lat
  i=1L
  while not eof(lun) do begin
    readu, inlun, time, n
    lon=fltarr(n)
    lat=fltarr(n)
    readu, inlun, lon, lat
    writeu, outlun, time, long(n), lon, lat
    i=i+1
  endwhile
  
  free_lun, inlun, outlun
  
  openu, outlun, newfile
  readu, outlun, magic
  writeu, outlun, i
  free_lun, outlun  
    
end
