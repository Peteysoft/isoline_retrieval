retbase="ecmwf1999030218_cld36s0.001SQ.ret.b"
truthfile="ecmwf2000101600_36s0.001SQ.cls"

openr, 1, retbase+".cls"
n=(fstat(1)).size/4
cls=lonarr(n)
readu, 1, cls
close, 1

openr, 1, retbase+".con"
con=fltarr(n)
readu, 1, con
close, 1

openr, 1, truthfile
truth=lonarr(n)
readu, 1, truth
close, 1

;divide into 0 and 1 results:
ind0=where(cls eq 0, n0)
ind1=where(cls eq 1, n1)

;convert confidence into cond. prob. for 1:
con0=0.5+con[ind0]/2
con1=0.5+con[ind1]/2

;normalize the classes:
clsn0=(truth[ind0] eq 0)/con0
clsn1=truth[ind1]/con1

sind0=sort(con0)
sind1=sort(con1)

clst0=total(clsn0[sind0], /cumulative)
clst1=total(clsn1[sind1], /cumulative)

;create the "confidence" axis:
ntick=6
contick=0.5*findgen(ntick)/(ntick-1)+0.5
xtick0=interpol(clst0, con0[sind0], contick)
xtick0[ntick-1]=clst0[n0-1]
xtick1=interpol(clst1, con1[sind1], contick)
xtick1[ntick-1]=clst1[n1-1]
ytick0=interpol(findgen(n0), con0[sind0], contick)
ytick0[ntick-1]=n0-1
ytick1=interpol(findgen(n1), con1[sind1], contick)
ytick1[ntick-1]=n1-1
tname=string(contick, format="(f4.1)")

set_plot, "ps"
device, /portrait, xsize=6.5, ysize=9, xoffset=1, yoffset=1, /inches

!p.multi=[0, 1, 2]

;window, 0
plot, clst0, findgen(n0), psym=3, xstyle=8, ystyle=8, xmargin=[7, 7]
oplot, findgen(10)*29000, findgen(10)*29000
axis, 0, !y.crange[1], /xaxis, xtickv=xtick0, xtickname=tname, xticks=ntick
axis, !x.crange[1], 0, /yaxis, ytickv=ytick0, ytickname=tname, yticks=ntick
;window, 1
plot, clst1, findgen(n1), psym=3, xstyle=8, ystyle=8, xmargin=[7, 7]
oplot, findgen(10)*29000, findgen(10)*29000
axis, 0, !y.crange[1], /xaxis, xtickv=xtick1, xtickname=tname, xticks=ntick
axis, !x.crange[1], 0, /yaxis, ytickv=ytick1, ytickname=tname, yticks=ntick

!p.multi=[0, 1, 1]
device, /close_file
device, /close

set_plot, "x"

end

