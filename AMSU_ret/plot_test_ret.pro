;build up the data structure:
t={date_rec, year: 0, pad: 0, day:0.D}

np=0L
nchan=0L
nscan=0L

openr, 1, "testret.dat"
comment=bytarr(500)
readu, 1, comment
t0={date_rec}
tf={date_rec}
readu, 1, t0, tf
readu, 1, nscan, np, nchan

rec={date:{date_rec}, lon: fltarr(np), lat:fltarr(np), bt:fltarr(nchan, np)}

data=replicate(rec, nscan)

readu, 1, data

close, 1

end
	
