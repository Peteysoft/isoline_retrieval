@rawread_amsu_l1c

;path1="../data/cret/20030101/"
path2="../data/cret/20020802/"

;spawn, "ls "+path1, filelist1, count=nf1
spawn, "ls "+path2, filelist2, count=nf2

plot, findgen(360)-180, findgen(360)/2-90, /nodata

;ind1=where(strmid(filelist1, 6, 1) eq 'M', cnt1)
;ind2=where(strmid(filelist2, 6, 1) eq 'M', cnt2)

r=findgen(21)*10
g=fltarr(21)
b=250-findgen(21)*10

cind=rgbtoindex(r, g, b)

cnt1=0
for i=0L, cnt1-1 do begin
  data=rawread_amsu_l1c(path1+filelist1[ind1[i]])
  ind=where(data.bt gt 0)
  oplot, (data.lon)[ind], (data.lat)[ind], psym=3
endfor

for i=0L, nf2-1 do begin
  data=rawread_amsu_l1c(path2+filelist2[i], missing=missing)
  ind=where(data.bt gt missing)
  ;round_sym, 0.1, /fill
  plot, findgen(360)-180, findgen(360)/2-90, /nodata
  color=fix((data.bt)[ind]*10+10)+1
  plots, (data.lon)[ind], (data.lat)[ind], psym=2, color=cind[color], symsize=0.002
  stop
endfor

end

