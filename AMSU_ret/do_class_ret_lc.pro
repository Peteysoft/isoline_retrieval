@time_lib

;start the collocation process:
tstart=convert_date('2002/12/31')
;tend=convert_date('2004/1/1')
;tend=convert_date('2002/1/1')

oneday={time_str, 0, 0, 1, 0, 0, 0}
;ndays=n_intervals(tstart, tend, oneday)+1
ndays=7
days=time_array(tstart, oneday, ndays)

basepath="/freax/storage/home/pmills/isoline_ret2/"

amsu_base="/talax/storage/noaa/"
satid=15

level3_afile="a.l1c"
corrected_afile_land="ac_land.l1c"
corrected_afile_sea="ac_sea.l1c"

level3_bfile="b.l1c"
corrected_bfile="bc.l1c"

sid_str=string(satid, format="(i2.2)")
acoeffile_land=basepath+"data/limb_corr/n"+sid_str+"limb_amsua_land.txt"
acoeffile_sea=basepath+"data/limb_corr/n"+sid_str+"limb_amsua_sea.txt"

if (satid eq 17) then begin
  bcoeffile=basepath+"data/limb_corr/n17limb_amsub_global.txt"
endif else begin
  bcoeffile=basepath+"data/limb_corr/n"+sid_str+"limb_amsub_land.txt"
endelse

border_base_land="iso0.00z36s0.001sq_land"
border_base_sea="iso0.00z36s0.001sq_sea"

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
    resultfile=basepath+"data/cret_lc/"+string(days[i].year, days[i].month, days[i].day, $
		format="(i4.4, i2.2, i2.2)")+"/CRET."+strmid(afiles[j], 9, 33)+".r1c"

    ;convert the level 1b files to level 1c:
    ;shell_command =  'ssh marvin "' + script + " " + $
    ;            amsu_base+amsua_path+afiles[j] +" " + level3_afile + '"'
    shell_command =  script + " " + amsu_base+amsua_path+afiles[j] +" " + level3_afile 
    spawn, shell_command, exit_status=err
    if err ne 0 then continue
    ;shell_command =  'ssh marvin "' + script + " " + $
    ;            amsu_base+amsub_path+bfiles[j] +" " + level3_bfile + '"'
    shell_command =  script + " " + amsu_base+amsub_path+bfiles[j] +" " + level3_bfile
    spawn, shell_command, exit_status=err
    if err ne 0 then continue
    ;stop

    ;apply limb corrections:
    spawn, "../src/amsu_lc "+level3_afile+" "+acoeffile_land+" "+corrected_afile_land
    spawn, "../src/amsu_lc "+level3_afile+" "+acoeffile_sea+" "+corrected_afile_sea
    spawn, "../src/amsu_lc "+level3_bfile+" "+bcoeffile+" "+corrected_bfile

    command="../src/retrieve_isoline_lc "+border_base_land+" "+border_base_sea+" "+ $
			corrected_afile_land+" "+corrected_afile_sea+" "+corrected_bfile + $
			" "+resultfile
    print, command
    spawn, command, exit_status=err
    if (err ne 0) then stop

  endfor
endfor

end

