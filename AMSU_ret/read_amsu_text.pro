function read_amsu_text, filename, interp_result=interp_result
  max_line=1000000

  header=""
  openr, lun, filename, /get_lun
  readf, lun, header
  readf, lun, header
  instr=strmid(header, 15, 1)
  if instr eq 'B' then begin
    nchan=5
  endif else if instr eq 'A' then begin
    nchan=15
  endif else begin
    print, "Routine only define for A and B instruments..."
    return, 0
  endelse
  readf, lun, header
  readf, lun, header

  if keyword_set(interp_result) then begin
    readf, lun, header
    readf, lun, header
    readf, lun, header
    readf, lun, header
  endif

  record={date:"", lon:0., lat:0., bt:fltarr(nchan)}
  data=replicate(record, max_line)

  on_ioerror, finish
  readf, lun, data, format="("+string(max_line)+"(a25, f0, f0, "+string(nchan)+"f0, /))"

  finish:
    ind=where(data.bt[0] ne 0 or data.date ne "")
    

  return, data[ind]

end
