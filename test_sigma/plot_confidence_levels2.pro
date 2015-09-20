pro plot_confidence_levels, nlev

  set_plot, "ps"
  device, filename="confidence_legend.ps"
  device, /landscape, xsize=7, ysize=1.5, xoffset=1, yoffset=10, /inches
  device, /color, bits_per_pixel=8

  ;create the color index:
  ncol=nlev
  ;green=[0, fltarr(ncol), 250]
  ;red=[0, findgen(ncol+1)*200/(ncol), 250]
  ;blue=[0, 200-findgen(ncol)*200/(ncol-1), 250]
  
  ;red=round([0, findgen(nlev)/(nlev-1)*155+50, 255, 200, 200])
  ;blue=round([0, 255-findgen(nlev)/(nlev-1)*255, 255, 0, 221])
  ;green=round([0, 150-findgen(nlev)/(nlev-1)*100, 255, 0, 50])

  red=findgen(nlev)*255/(nlev-1)
  blue=255-red
  green=fltarr(nlev)
  
  tvlct, red, green, blue, 1
  colind=rgbtoindex(red, green, blue)
  colind=indgen(nlev+2)
  ccol=colind[0:nlev-1]

  conlev=findgen(nlev)/(nlev-1)

  legend=rebin(conlev, nlev, 2)

  ;window, 0, ysize=500, xsize=150

  contour, legend, conlev, [0, 1], levels=conlev, xtitle="Confidence", $
	  charsize=1.3, c_colors=colind[1:nlev], /fill, $
	charthick=2, ystyle=4;, $
;	  ymargin=[2.5, 2.5], xmargin=[7, 2]
	  ;background=colind[nlev+1], color=colind[0]

  device, /close
  set_plot, "x"

  ;image=tvrd(0)
  ;ind=where(image eq 250)
  ;image=fix(round(image*(nlev-1.)/200))
  ;image[ind]=nlev+1

  ;write_tiff, "confidence_legend.tif", image, $
;	  	red=fix(red), blue=fix(blue), green=fix(green)

  stop

end
	      
