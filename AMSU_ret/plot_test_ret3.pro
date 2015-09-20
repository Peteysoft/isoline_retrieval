;path1="../data/cret/20030101/"
path2="../data/cret/20030102/"

;spawn, "ls "+path1, filelist1, count=nf1
spawn, "ls "+path2, filelist2, count=nf2

set_plot, "ps"
device, /landscape, xsize=9.5, ysize=6.5, xoffset=1, yoffset=10.5, /inches
device, filename="testret3.ps"
device, /color, bits_per_pixel=8

plot, findgen(360)-180, findgen(360)/2-90, /nodata

;ind1=where(strmid(filelist1, 6, 1) eq 'M', cnt1)
;ind2=where(strmid(filelist2, 6, 1) eq 'M', cnt2)

r=[0, 250, 0, 0, 255]
g=[0, 0, 250, 0, 255]
b=[0, 0, 0, 250, 255]

tvlct, r, g, b

red=1
green=2
blue=3

cnt1=0
for i=0L, cnt1-1 do begin
  data=rawread_amsu_l1c(path1+filelist1[ind1[i]])
  ind=where(data.bt gt 0)
  oplot, (data.lon)[ind], (data.lat)[ind], psym=3
endfor

for i=0L, nf2-1 do begin
  data=rawread_amsu_l1c(path2+filelist2[i])
  ind=where(data.bt gt 0)
  sid=strmid(filelist2[i], 6, 1)
  if (sid eq 'K') then color=red else if sid eq 'L' then color=green else color=blue
  oplot, (data.lon)[ind], (data.lat)[ind], psym=3, color=color
endfor

device, /close_file
device, /close
set_plot, "x"

end

