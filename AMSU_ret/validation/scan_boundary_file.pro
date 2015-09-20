@time_lib

function scan_boundary_file, filename

  openr, lun, filename, /get_lun

  magic=0L
  nrec=0L
  readu, lun, magic, nrec
  print, nrec, " records found in file, ", filename

  date=replicate({time_str}, nrec)
  n=0L
  date1={time_str}

  for i=0L, nrec-1 do begin

    readu, lun, date1, n
    date[i]=date1

    print, n, " ", time_string(date1)

    point_lun, -lun, cur
    point_lun, lun, cur+8*n
  endfor

  stop
  free_lun, lun

  return, date

end

