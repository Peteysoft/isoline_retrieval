function rawread_amsu_l1c, filename, missing=missing

  ;build up the data structure:
  t={date_rec, year: 0, pad: 0, day:0.D}

  np=0L
  nchan=0L
  nscan=0L
  missing=0.

  openr, lun, filename, /get_lun
  comment=bytarr(500)
  readu, lun, comment
  print, string(comment)
  t0={date_rec}
  tf={date_rec}
  readu, lun, t0, tf
  readu, lun, nscan, np, nchan, missing

  rec={date:{date_rec}, lon: fltarr(np), lat:fltarr(np), bt:fltarr(nchan, np)}

  data=replicate(rec, nscan)

  readu, lun, data

  free_lun, lun

  return, data

end
	
