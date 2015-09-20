
#include "traj_int_obj.h"
#include "colloc_obj.h"

class tcolloc_obj: public colloc_obj {

  protected:
    //for integrating trajectories:
    traj_int_obj *tion;		//N. hemisphere
    traj_int_obj *tios;		//S. hemisphere

    interpol_index dt;		//time step
    long nt;			//number of steps

  public:
    //for passing back the trajectory coordinates:
    float traj_lon, traj_lat;
 
    //amsu_tinterp_obj();
    tcolloc_obj(char *indexfile, long mnl, char *basepath, char *nfile, char *sfile, 
		interpol_index dt1, long nt1);
    ~tcolloc_obj();

    void set_page_size(long page_size);

    virtual float colloc_forward(float lon,               //location of test point
                float lat,
                time_class t0,          //time of test point
                float critdist,
                time_class &tforward,
                long &pixnum,            //pixel number
                char *&filename,         //(pointer to) filename
                long &nchan,            //total number of channels
                float *&val);            //(pointer to) values of nearest neighbour

    virtual float colloc_backward(float lon,              //location of test point
                float lat,
                time_class t0,          //time of test point
                float critdist,
                time_class &tbackward,
                long &pixnum,            //pixel number
                char *&filename,         //(pointer to) filename
                long &nchan,            //total number of channels
                float *&val);            //(pointer to) values of nearest neighbour

};

