

;start the collocation process:
tstart=convert_date('2002/12/31')
;tend=convert_date('2004/1/1')
;tend=convert_date('2002/1/1')

oneday={time_str, 0, 0, 1, 0, 0, 0}
;ndays=n_intervals(tstart, tend, oneday)+1
ndays=7
days=time_array(tstart, oneday, ndays)

basepath="/freax/storage/home/pmills/isoline_ret2/data/"

for i=0L, ndays-1 do begin
  amsua_path=string("amsul1c/amsua", days[i].year, days[i].month, days[i].day, $
   		format="(a, i4.4, i2.2, i2.2)")+"/"
  amsub_path=string("amsul1c/amsub", days[i].year, days[i].month, days[i].day, $
   		format="(a, i4.4, i2.2, i2.2)")+"/"
  spawn, "ls "+basepath+amsua_path, afiles, count=na
  spawn, "ls "+basepath+amsub_path, bfiles, count=nb

  for j=0L, na-1 do begin
    blen=strlen(afiles[j])
    check_bfile=afiles[j]
    strput, check_bfile, 'B', 6
    if (bfiles[j] ne check_bfile) then begin
      print, "A and B names do not agree: ", afiles[j], bfiles[j]
      stop
    endif
    resultfile=basepath+"cret/"+string(days[i].year, days[i].month, days[i].day, $
		format="(i4.4, i2.2, i2.2)")+"/CRET."+strmid(afiles[j], 9, 33)+".r1c"
    ;stop
    fin1=basepath+amsua_path+afiles[j]
    fin2=basepath+amsub_path+bfiles[j]
    command="../src/retrieve_isoline2 "+fin1+" "+fin2+" "+resultfile
    print, command
    spawn, command

  endfor
endfor

end

