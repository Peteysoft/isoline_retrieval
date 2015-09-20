;@traj_utils

function distance, p1, p2

;+*************************************************************************
;		DISTANCE
;**************************************************************************
;
; usage:	d=distance(v1, v2)
;
; purpose:	Returns the approximate distance (correct in the limit as 
;		it approaches 0) between two points defined in a spherical 
;		polar coords with constant radius.  For the poles, converts 
;		to a locally Cartesian system.
;
; parameters:	v1		The first vector.  Should be a two element vector
;				containing longitude and latitude respectively.
;
;		v2		The second vector in the same format as the first.
;
; history:	Formally documented, 2002-02-14 Peter Mills (peter.mills@nrl.navy.mil)
;
;-**************************************************************************

  mperdeg=40000000./360
  latthresh=87.5

  if abs(p1[1]) ge latthresh or abs(p2[1]) ge latthresh then begin
  ;if points are very close to the pole, use a local cartesian
  ;system calculate the distance:
    rearth=20000000/!pi
    r1=rearth*!pi*(90-abs(p1[1]))/180
    r2=rearth*!pi*(90-abs(p2[1]))/180
    x1=r1*cos(!pi*p1[0]/180)
    y1=r1*sin(!pi*p1[0]/180)
    x2=r2*cos(!pi*p2[0]/180)
    y2=r2*sin(!pi*p2[0]/180)
    return, sqrt((x1-x2)^2+(y1-y2)^2)
  endif else begin
    dlon=abs(p1[0]-p2[0])
    if dlon gt 180 then dlon=360-dlon
;    print, dlon
    dx=dlon*cos(!pi*(p1[1]+p2[1])/360)
    dy=p1[1]-p2[1]

    return, mperdeg*sqrt(dx^2+dy^2)
  endelse
end

function midpoint, p1, p2

;+*************************************************************************
;		MIDPOINT
;**************************************************************************
;
; usage:	v=midpoint(v1, v2)
;
; purpose:	Returns the midpoint between two points specified on a 
;		spherical-polar coordinate system with constant radius.  
;		Correct in the limit as the distance approaches zero.  For 
;		the poles, converts to a locally Cartesian system.
;
; parameters:	v1		The first vector.  Should be a two element vector
;				containing longitude and latitude respectively.
;
;		v2		The second vector in the same format as the first.
;
; history:	Formally documented, 2002-02-14 Peter Mills (peter.mills@nrl.navy.mil)
;
;-**************************************************************************


  if abs(p1[1]) ge 87.5 or abs(p2[1]) ge 87.5 then begin
  ;if points are very close to the pole, use a local cartesian
  ;system to do the bisection:
    r1=!pi*(90-abs(p1[1]))/180
    r2=!pi*(90-abs(p2[1]))/180
    x1=r1*cos(!pi*p1[0]/180)
    y1=r1*sin(!pi*p1[0]/180)
    x2=r2*cos(!pi*p2[0]/180)
    y2=r2*sin(!pi*p2[0]/180)
    xnew=(x1+x2)/2
    ynew=(y1+y2)/2
    lon=180*atan(ynew, xnew)/!pi
    if lon lt 0 then lon=360+lon
    lat=abs(p1[1])/p1[1]*(90-180*sqrt(xnew^2+ynew^2)/!pi)
;    stop
    return, [lon, lat]
  endif else begin
  
    if abs(p1[0]-p2[0]) gt 180 then lonnew=(p1[0]+p2[0]-360)/2 else lonnew=(p1[0]+p2[0])/2
    if lonnew lt 0 then lonnew=360+lonnew
;    if lonnew gt 360. then lonnew=lonnew-360.
;    stop

    return, [lonnew, (p1[1]+p2[1])/2]
  endelse
end

function point_spacing, x, y, nowrap=nowrap
;+*************************************************************************
;		POINT_SPACING
;**************************************************************************
;
; usage:	d=point_spacing(x, y [, nowrap])
;
; purpose:	Given a set of coordinate pairs which are assumed to define
;		a closed contour or mixing boundary, returns the distances
;		between adjacent points, include first and last.
;
; parameters:	x		An array of values containing the first coordinate.
;
;		x		An array of values containing the second coordinate.
;
; keywords:	nowrap	Doesn't return the distance between the first and 
;				last points.
;
; dependencies:		distance
;
; history:	Formally documented, 2002-02-14 Peter Mills (peter.mills@nrl.navy.mil)
;		-2002-3-21 PM:  added nowrap keyword
;
;-**************************************************************************

  n=n_elements(x)
  if keyword_set(nowrap) then distance=fltarr(n-1) else distance=fltarr(n)
  for i=0L, n-2 do begin
    distance[i]=distance([x[i], y[i]], [x[i+1], y[i+1]])
  endfor
  if not keyword_set(nowrap) then distance[n-1]=distance([x[n-1], y[n-1]], [x[0], y[0]])
  return, distance
end


pro generate_circle, centre, radius, n, x, y
;+*************************************************************************
;		GENERATE_CIRCLE
;**************************************************************************
;
; usage:	generate_circle, centre, radius, n, x, y
;
; purpose:	Generates an approximate circle on the surface of a sphere.
;		Not for use near the poles.
;
; parameters:	centre	A two element vector containing the longitude and
;				latitude of the centre of the circle.
;
;		radius	The radius of the circle in metres.
;
;		n		The number of points in the circle.
;
;		x		An array of values containing the longitudes of the
;				circle.
;
;		x		An array of values containing the latitudes of the
;				circle.
;
; history:	Formally documented, 2002-02-14 Peter Mills (peter.mills@nrl.navy.mil)
;
;-**************************************************************************

  mperdeg=40000000./360
  angle=2*!pi*findgen(n)/n
  y=radius*sin(angle)/mperdeg+centre[1]
  x=radius*cos(angle)/mperdeg/cos(!pi*(y+centre[1])/360)+centre[0]
  dum=x
  
end

