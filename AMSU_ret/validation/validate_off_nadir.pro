@time_lib

;this program validates RTTOV simulations for off nadir case

RDRYAIR=287
RWATERVAPOUR=461
GACC=9.8

zangle=8800./89-50

satid='16'

n=200L

;a bunch of random dates (no repeats...)
t0=convert_date('2002/1/1')
oneday=convert_date('0/0/1')
nt=365*2

tarr=time_array(t0, oneday, nt)
sind=sort(randomu(seed, nt))
t=tarr[sind]

level3_afile="a.l1c"
level3_bfile="b.l1c"

script="/freax/storage/home/pmills/atovs_tools/shell_scripts/zamsu2l1c.sh"

field=["PP", "TT", "SQ", "CL", "CI", "CC", "UU", "VV"]

nlev=60

nlon=240
nlat=121

mask=intarr(240, 121)
openr, 1, "landmask.txt"
readf, 1, mask
close, 1

bta_nadir=fltarr(15, n)
btb_nadir=fltarr(5, n)
bta_off=fltarr(15, n)
btb_off=fltarr(5, n)

openw, 1, "random_off_nadir_direct.txt"
openw, 2, "random_nadir_direct.txt"

k=0
for i=0L, n-1 do begin
  print, i
  yrstr=string(t[j].year, format="(i4.4)")
  
  patha="/talax/storage/noaa/noaa"+satid+"_amsua_"+yrstr+ $
		string(t[j].month, t[j].day, format="('/', i2.2, '/', i2.2)")
  pathb="/talax/storage/noaa/noaa"+satid+"_amsub_"+yrstr+ $
		string(t[j].month, t[j].day, format="('/', i2.2, '/', i2.2)")

  spawn, "ls "+patha, flista
  spawn, "ls "+pathb, flistb

  nf=n_elements(flista)

  wf=fix(randomu(seed, 1)*nf)

  ;pick a random file for the date:
  fa=flista[wf]
  fb=flistb[wf]
  fb1=fb
  strput, fb1, 'A', 6
  if strmid(fa, 0, 42) ne strmid(fb1, 0, 42) then begin
    fb=""
    for j=0, n_elements(flistb)-1 do begin
      fb1=flistb[j]
      strput, fb1, 'A', 6
      if strmid(fa, 0, 42) eq strmid(fb1, 0, 42) then begin
        fb=flistb[j]
        break
      endif
    endfor
    if fb eq "" then stop
    i=i-1
    continue
  endif

  ;convert the data:
  ;convert the level 1b files to level 1c:
  ;shell_command =  'ssh marvin "' + script + " " + $
  ;		patha+"/"+fa +" " + level3_afile + '"'
  shell_command =  script + " " + patha+"/"+fa +" " + level3_afile
  spawn, shell_command, exit_status=err
  if err ne 0 then begin
    print, "Failed to convert file: ", fa
    stop
    i=i-1
    continue
  endif
  ;shell_command =  'ssh marvin "' + script + " " + $
  ;            pathb+"/"+fb +" " + level3_bfile + '"'
  shell_command =  script + " " + pathb+"/"+fb +" " + level3_bfile
  spawn, shell_command, exit_status=err
  if err ne 0 then begin
    print, "Failed to convert file: ", fb
    stop
    i=i-1
    continue
  endif
  ;spawn, "../../src/fswap_endian "+level3_afile
  ;spawn, "../../src/fswap_endian "+level3_bfile

  ;read in the data:
  adata=read_amsu_l1c(level3_afile)
  bdata=read_amsu_l1c(level3_bfile)

  ;select a random scan line:
  sla=randomu(seed, 1)*n_elements(adata)
  date=adata[sla].date
  slb=bin_search(bdata.date, adata[sla].date, /nearest)

  ;pick a random side of the scan line:
  lr=fix(randomu(seed, 1)*2)
  if lr eq 1 then begin
    lon=adata[sla].lon[29]
    lat=adata[sla].lat[29]
    bta_off[*, i]=adata[sla].bt[*, 29]
    btb_off[*, i]=bdata[slb].bt[*, 88]
  endif else begin
    lon=adata[sla].lon[0]
    lat=adata[sla].lat[0]
    bta_off[*, i]=adata[sla].bt[*, 0]
    btb_off[*, i]=bdata[slb].bt[*, 1]
  endelse

  ;get the appropriate ecmwf profile:
  if lon lt 0 then lon1=360+lon else lon1=lon
  get_ecmwf_profile, field, lon1, lat, date, profile
  profile[*, 0]=profile[*, 0]*100

  ;calculate thickness of first layer:
  wvp=profile[nlev-1, 2]*profile[nlev-1, 0]*1.60627
  req1=profile[nlev-1, 0]/((profile[nlev-1, 0]-wvp)/RDRYAIR+wvp/RWATERVAPOUR)

  wvp=profile[nlev-2, 2]*profile[nlev-2, 0]*1.60627
  req2=profile[nlev-2, 0]/((profile[nlev-2, 0]-wvp)/RDRYAIR+wvp/RWATERVAPOUR)

  c=(req1*profile[nlev-1, 1]+req2*profile[nlev-2, 1])/GACC/2
  z1=c*alog(profile[nlev-1, 0]/profile[nlev-2, 0])

  frac=2/z1

  p2m=profile[nlev-1, 0]*(1-frac)+profile[nlev-2, 0]*frac
  t2m=profile[nlev-1, 1]*(1-frac)+profile[nlev-2, 1]*frac
  q2m=profile[nlev-1, 2]*(1-frac)+profile[nlev-2, 2]*frac

  u2m=profile[nlev-1, 6]*(1-frac)+profile[nlev-2, 6]*frac
  v2m=profile[nlev-1, 7]*(1-frac)+profile[nlev-2, 7]*frac

  ;get land/sea flag:
  lonind=lon1/1.5
  latind=(90-lat)/1.5
  ls=interpolate(mask, lonind, latind) gt 0.5

  printf, 1, lon, lat, ls, profile[nlev-1, 1], profile[nlev-1, 0], $
		t2m, q2m, u2m, v2m, format="(9g16.6)"
  printf, 1, profile[*, 1:5], format="(30(10g16.6/), $)"
  printf, 1, fltarr(nlev, 2), format="(12(10g16.6/), $)"

  ;stop

  ;same thing for nadir case:
  ;pick a random side of the scan line:
  lr=fix(randomu(seed, 1)*2)
  if lr eq 1 then begin
    lon=adata[sla].lon[15]
    lat=adata[sla].lat[15]
    bta_nadir[*, i]=adata[sla].bt[*, 15]
    btb_nadir[*, i]=bdata[slb].bt[*, 45]
  endif else begin
    lon=adata[sla].lon[14]
    lat=adata[sla].lat[14]
    bta_nadir[*, i]=adata[sla].bt[*, 14]
    btb_nadir[*, i]=bdata[slb].bt[*, 44]
  endelse

  ;get the appropriate ecmwf profile:
  if lon lt 0 then lon1=360+lon else lon1=lon
  get_ecmwf_profile, field, lon1, lat, date, profile
  profile[*, 0]=profile[*, 0]*100

  ;calculate thickness of first layer:
  wvp=profile[nlev-1, 2]*profile[nlev-1, 0]*1.60627
  req1=profile[nlev-1, 0]/((profile[nlev-1, 0]-wvp)/RDRYAIR+wvp/RWATERVAPOUR)

  wvp=profile[nlev-2, 2]*profile[nlev-2, 0]*1.60627
  req2=profile[nlev-2, 0]/((profile[nlev-2, 0]-wvp)/RDRYAIR+wvp/RWATERVAPOUR)

  c=(req1*profile[nlev-1, 1]+req2*profile[nlev-2, 1])/GACC/2
  z1=c*alog(profile[nlev-1, 0]/profile[nlev-2, 0])

  frac=2/z1

  p2m=profile[nlev-1, 0]*(1-frac)+profile[nlev-2, 0]*frac
  t2m=profile[nlev-1, 1]*(1-frac)+profile[nlev-2, 1]*frac
  q2m=profile[nlev-1, 2]*(1-frac)+profile[nlev-2, 2]*frac

  u2m=profile[nlev-1, 6]*(1-frac)+profile[nlev-2, 6]*frac
  v2m=profile[nlev-1, 7]*(1-frac)+profile[nlev-2, 7]*frac

  ;get land/sea flag:
  lonind=lon1/1.5
  latind=(90-lat)/1.5
  ls=interpolate(mask, lonind, latind) gt 0.5

  printf, 2, lon, lat, ls, profile[nlev-1, 1], profile[nlev-1, 0], $
		t2m, q2m, u2m, v2m, format="(9g16.6)"
  printf, 2, profile[*, 1:5], format="(30(10g16.6/), $)"
  printf, 2, fltarr(nlev, 2), format="(12(10g16.6/), $)"

endfor

close, 1
close, 2

openw, 1, "btA_off.txt"
printf, 1, bta_off
close, 1

openw, 1, "btB_off.txt"
printf, 1, btb_off
close, 1

openw, 1, "btA_nadir.txt"
printf, 1, bta_nadir
close, 1

openw, 1, "btB_nadir.txt"
printf, 1, btb_nadir
close, 1

end

