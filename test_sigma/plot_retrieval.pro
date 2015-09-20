@ecmwf_iso_plotting

date='1999030218'

sqthresh=0.001
sigma=36

clsfile="ecmwf"+date+"_cld36s0.001SQ0.00z.ret.b.cls"
confile="ecmwf"+date+"_cld36s0.001SQ0.00z.ret.b.con"

nlon=240
nlat=121

cls=lonarr(nlon, nlat)
openr, lun, clsfile, /get_lun
readu, lun, cls
free_lun, lun

con=fltarr(nlon, nlat)
openr, lun, confile, /get_lun
readu, lun, con
free_lun, lun

ecmwf_iso_plot_comp, '1999/3/2-18', cls, sigma, sqthresh, con, $
		psfile='ret'+date+'.ps'
ecmwf_iso_plot_range, '1999/3/2-18', cls, sigma, sqthresh, con, 0.81, $
		psfile='con'+date+'.ps'

end
