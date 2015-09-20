
;read the confidence ratings at all the training points:
datelist=['1998082718', '1999030218', '1999092212', '2000101600', $
		'2001122006', '2002061212']
flist="ecmwf"+datelist+"_cld36s0.001SQ0.00z.ret.b.con"
nf=n_elements(flist)
npt=240L*121
con=fltarr(npt*nf)

confield=fltarr(npt)

for i=0, nf-1 do begin
  openr, 1, flist[i]
  readu, 1, confield
  close, 1
  con[i*npt]=confield
endfor

;transform to 1-P(c|x), where c is the "winning" class:
r=1-(con+1)/2

nd=6

sind=sort(con)

set_plot, "ps"
device, filename="confunc_the2.ps"
device, /portrait, xsize=6.5, ysize=6.5, /inches

delta=total(r[sind], /cum)

plot, con[sind], delta/delta[npt*nf-1], title="Tolerance vs. confidence: theoretical", $
		xtitle="confidence rating", ytitle="tolerance interval"

for i=1, nd-1 do begin
  alpha=1+1.*i/nd

  delta=total(r[sind]^(1/alpha), /cum)

  oplot, con[sind], delta/delta[npt*nf-1]

endfor

device, /close
device, /close_file

set_plot, "x"

end

