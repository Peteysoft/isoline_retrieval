
;read the confidence ratings at all the training points:
n=85965L
con=fltarr(n)

openr, 1, "p2ndorder.con"
readu, 1, con
close, 1

;transform to 1-P(c|x), where c is the "winning" class:
r=(1-con)/2

nd=6

sind=sort(con)

set_plot, "ps"
device, filename="confunc_the.ps"
device, /portrait, xsize=6.5, ysize=6.5, /inches

delta=total(r[sind], /cum)

plot, con[sind], delta/delta[n-1], title="Tolerance vs. confidence: theoretical", $
		xtitle="confidence rating", ytitle="tolerance interval"
ind=interpol(findgen(n), con[sind], 0.5)
xyouts, con[sind[ind]], delta[ind]/delta[n-1], string(1., format="(f3.1)"), /data

for i=1, nd-1 do begin
  alpha=1+1.*i/(nd-1)

  delta=total(r[sind]^(1/alpha), /cum)

  oplot, con[sind], delta/delta[n-1]

  ind=interpol(findgen(n), con[sind], 0.5)
  xyouts, con[sind[ind]], delta[ind]/delta[n-1], string(alpha, format="(f3.1)"), /data

endfor

device, /close
device, /close_file

set_plot, "x"

end

