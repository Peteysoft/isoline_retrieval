#include <math.h>

#include "amsu_interp_base.h"

#define NORMAL 0
#define RAW 1

#define BIN_NLON 360
#define BIN_NLAT 180

struct amsu_index_bin {
  long n;
  time_class *date;
  long *file_index;
  long *scan_index;
};

class amsu_interp_obj: public amsu_interp_base {

  protected:
    //the index:
    long *bin_nlon, bin_nlat;
    float polar_latthresh;
    amsu_index_bin **bin;

    //calculate bin indices:
    int binind(float lon, float lat, long &lonind, long &latind);

    float lon_offset;

  public:
    amsu_interp_obj();
    amsu_interp_obj(char *initfile, long mnl, char *basepath=NULL);
    ~amsu_interp_obj();

    virtual int int_forward(float lon,		//location of test point
		float lat,
		time_class t0, 		//time of test point
		float critdist,
		time_class &tforward,
		long &nchan,		//total number of channels
		float *val[4],		//values at eight interpolation neighbours
		float weight[4],	//their weights
		int no_extrap=1);

    virtual int int_backward(float lon,		//location of test point
		float lat,
		time_class t0, 		//time of test point
		float critdist,
		time_class &tbackward,
		long &nchan,		//total number of channels
		float *val[4],		//values at eight interpolation neighbours
		float weight[4],	//their weights
		int no_extrap=1);

    virtual float colloc_forward(float lon,		//location of test point
		float lat,
		time_class t0, 		//time of test point
		float critdist,
		time_class &tforward,
		long &scanline,		//scan number
		long &pixnum,		//pixel number
		char *&filename,		//(pointer to) filename
		long &nchan,		//total number of channels
		float *&val);		//(pointer to) values of nearest neighbour

    virtual float colloc_backward(float lon,		//location of test point
		float lat,
		time_class t0, 		//time of test point
		float critdist,
		time_class &tbackward,
		long &scanline,		//scan number
		long &pixnum,		//pixel number
		char *&filename,		//(pointer to) filename
		long &nchan,		//total number of channels
		float *&val);		//(pointer to) values of nearest neighbour

};

