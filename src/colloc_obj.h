#include <math.h>

#include "colloc_base.h"

#define NORMAL 0
#define RAW 1

class colloc_obj: public colloc_base {

  protected:
    //indexes each entry by date:
    short *tind;
    int32_t nind;
    //time as a floating point to save space:
    time_class tref;
    float *t;

    int index_data();
    int read_flist(char *initfile);
    int read_indexfile(char *initfile);
    int write_indexfile(char *filename);

  public:
    colloc_obj();
    colloc_obj(char *initfile, long mnl, char *basepath1=NULL);
    ~colloc_obj();

    virtual float colloc_forward(float lon,		//location of test point
		float lat,
		time_class t0, 		//time of test point
		float critdist,
		time_class &tforward,
		long &pixnum,		//pixel number
		char *&filename,		//(pointer to) filename
		long &nchan,		//total number of channels
		float *&val);		//(pointer to) values of nearest neighbour

    virtual float colloc_backward(float lon,		//location of test point
		float lat,
		time_class t0, 		//time of test point
		float critdist,
		time_class &tbackward,
		long &pixnum,		//pixel number
		char *&filename,		//(pointer to) filename
		long &nchan,		//total number of channels
		float *&val);		//(pointer to) values of nearest neighbour

};

