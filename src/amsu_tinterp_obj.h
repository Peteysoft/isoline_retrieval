
#include "traj_int_obj.h"
#include "amsu_interp_obj.h"

class amsu_tinterp_obj: public amsu_interp_obj {

  protected:
    //for integrating trajectories:
    traj_int_obj *tion;		//N. hemisphere
    traj_int_obj *tios;		//S. hemisphere

    interpol_index dt;		//time step
    long nt;			//number of steps
 
  public:
    //amsu_tinterp_obj();
    amsu_tinterp_obj(char *indexfile, long mnl, char *basepath, char *nfile, char *sfile, 
		interpol_index dt1, long nt1);
    ~amsu_tinterp_obj();

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

    virtual float colloc_forward(float lon,               //location of test point
                float lat,
                time_class t0,          //time of test point
                float critdist,
                time_class &tforward,
                long &scanline,          //scan number
                long &pixnum,            //pixel number
                char *&filename,         //(pointer to) filename
                long &nchan,            //total number of channels
                float *&val);            //(pointer to) values of nearest neighbour

    virtual float colloc_backward(float lon,              //location of test point
                float lat,
                time_class t0,          //time of test point
                float critdist,
                time_class &tbackward,
                long &scanline,          //scan number
                long &pixnum,            //pixel number
                char *&filename,         //(pointer to) filename
                long &nchan,            //total number of channels
                float *&val);            //(pointer to) values of nearest neighbour

};

