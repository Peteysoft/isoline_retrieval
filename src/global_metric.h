
//approximate squared distance for lon-lat coordinates:
float sdist(float lon1, float lat1, float lon2, float lat2);

//interpolate in longitude and latitude:
void interpolate_lon_lat(float lon1, float lat1, float lon2, float lat2, 
		double frac, float &lon, float &lat);

