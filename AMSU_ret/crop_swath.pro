@time_lib

;start the collocation process:
tstart=convert_date('2002/12/31')
;tend=convert_date('2004/1/1')
;tend=convert_date('2002/1/1')

oneday={time_str, 0, 0, 1, 0, 0, 0}
;ndays=n_intervals(tstart, tend, oneday)+1
ndays=7
days=time_array(tstart, oneday, ndays)

ncrop=3

basepath="/freax/storage/home/pmills/isoline_ret2/"

for i=0L, ndays-1 do begin
  path1=basepath+"data/cret/"+string(days[i].year, days[i].month, days[i].day, $
		format="(i4.4, i2.2, i2.2)")+"/"
  path2=basepath+"data/cret_cr/"+string(days[i].year, days[i].month, days[i].day, $
		format="(i4.4, i2.2, i2.2)")+"/"

  spawn, "ls "+path1, file, count=nf

  for j=0L, nf-1 do begin
    command="../src/crop_swath -r "+path1+file[j]+" "+string(ncrop)+" "+path2+file[j]
    print, command
    spawn, command

  endfor
endfor

end

