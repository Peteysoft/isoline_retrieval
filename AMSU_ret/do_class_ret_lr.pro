@time_lib

;start the collocation process:
tstart=convert_date('2001/1/1')
;tend=convert_date('2004/1/1')
;tend=convert_date('2002/1/1')

oneday={time_str, 0, 0, 1, 0, 0, 0}
;ndays=n_intervals(tstart, tend, oneday)+1
ndays=7
days=time_array(tstart, oneday, ndays)

basepath="/freax/storage/home/pmills/isoline_ret2/"

amsu_base="/talax/storage/noaa/"
;amsu_base="/freax/storage/noaa/"
satid=16

level3_afile=basepath+"ret/a.l1c"
level3_bfile=basepath+"ret/b.l1c"

script="/freax/storage/home/pmills/atovs_tools/shell_scripts/zamsu2l1c.sh"

for i=0L, ndays-1 do begin
  amsua_path=string("noaa", satid, "_amsua_", days[i].year, "/", days[i].month, "/", days[i].day, $
   		format="(a, i2.2, a, i4.4, a1, i2.2, a1, i2.2)")+"/"
  amsub_path=string("noaa", satid, "_amsub_", days[i].year, "/", days[i].month, "/", days[i].day, $
   		format="(a, i2.2, a, i4.4, a1, i2.2, a1, i2.2)")+"/"
  spawn, "ls "+amsu_base+amsua_path, afiles, count=na
  spawn, "ls "+amsu_base+amsub_path, bfiles, count=nb

  for j=0L, na-1 do begin
    blen=strlen(afiles[j])
    check_bfile=afiles[j]
    strput, check_bfile, 'B', 6
    if (bfiles[j] ne check_bfile) then begin
      print, "A and B names do not agree: ", afiles[j], bfiles[j]
      stop
    endif
    resultfile=basepath+"data/cret_lr/"+string(days[i].year, days[i].month, days[i].day, $
		format="(i4.4, i2.2, i2.2)")+"/CRET."+strmid(afiles[j], 9, 33)+".r1c"

    ;convert the level 1b files to level 1c:
    shell_command =  'ssh marvin "' + script + " " + $
                amsu_base+amsua_path+afiles[j] +" " + level3_afile + '"'
    ;shell_command =  script + " " + amsu_base+amsua_path+afiles[j] +" " + level3_afile 
    spawn, shell_command, exit_status=err
    if err ne 0 then continue
    shell_command =  'ssh marvin "' + script + " " + $
                amsu_base+amsub_path+bfiles[j] +" " + level3_bfile + '"'
    ;shell_command =  script + " " + amsu_base+amsub_path+bfiles[j] +" " + level3_bfile
    spawn, shell_command, exit_status=err
    if err ne 0 then continue
    ;stop
    command="../src/retrieve_isoline_lr -s "+level3_afile+" "+level3_bfile+" "+resultfile
    print, command
    spawn, command

  endfor
endfor

end

