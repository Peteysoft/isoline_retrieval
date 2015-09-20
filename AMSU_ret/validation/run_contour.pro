;spawn, "ls c???.bev", filelist
spawn, "ls c_ecmwf???.bev", filelist

nwindfile="/home/pmills/isoline_ret2/data/vfield_2002_S36_N.ds"
swindfile="/home/pmills/isoline_ret2/data/vfield_2002_S36_S.ds"

t0='2002/9/3'
dt=0.1
nrk=3
n=200

initfile="cinit.cinit"
command="/home/pmills/bin/contour2 "+initfile

for i=0L, n_elements(filelist)-1 do begin
  outfile="smooth/"+strmid(filelist[i], 0, 10)+".result.bev"

  openw, 1, initfile
  printf, 1, nwindfile
  printf, 1, swindfile
  printf, 1, filelist[i]
  printf, 1, outfile
  printf, 1, t0
  printf, 1, dt
  printf, 1, n
  printf, 1, nrk
  printf, 1, 1.
  printf, 1, 10.
  close, 1

  spawn, command, exit_status=err
  if err ne 0 then stop

endfor

end
