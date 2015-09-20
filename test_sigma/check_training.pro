openr, 1, "train.vec"
nvar=0L
readu, 1, nvar
n=99965
data=fltarr(6, n)
readu, 1, data

k=0
for i=0L, n-1 do begin
  print, i
  for j=i+1, n-1 do begin
    if array_equal(data[*, i], data[*, j]) then begin
      print, k, ": vectors", i, " and", j, " are the same"
      k=k+1
    endif
  endfor
endfor

end

