#include <math.h>

#include "read_amsu.h"

#define NORMAL 0
#define RAW 1

#define BIN_NLON 360
#define BIN_NLAT 180

//once a point has been located, find three others and compute the weights:
void derive_amsu_int_weights(amsu_1c_data *data,//amsu data
                        long inds,              //the scan index
                        long indp,              //swath index
                        float lon,              //longitude of point
                        float lat,              //latitude of the point
                        float *val[4],          //values at each point
                        float weight[4]);       //returned weights

//base class (doesn't actually do anything...)
//contains those data and methods common to both flavours
class amsu_interp_base {

  protected:
    //list of files and data:
    long nfiles;		//total number of files
    long nloaded;		//number loaded
    long maxload;		//maximum number of loaded files allowed
    long offset;		//zero point of loaded data

    char **filelist;		//list of files in order of start date
    //time_class *tstart;	//start date of each file
    //time_class *tend;		//end date of each file
    amsu_1c_data **amsu_data;	//amsu data

    //file-swapping algorithm:
    amsu_1c_data *get_data_forward(long index);
    amsu_1c_data *get_data_backward(long index);

    //read routines:
    int swap_end_flag;
    long (*amsu_head_reader) (char *, time_class &, time_class &, long &, long &, int);
    amsu_1c_data * (*amsu_reader) (char *, int);

  public:
    virtual int int_forward(float lon,          //location of test point
		float lat,
		time_class t0,          //time of test point
		float critdist,
		time_class &tforward,
		long &nchan,            //total number of channels
		float *val[4],          //values at eight interpolation neighbours
		float weight[4],       //their weights
		int no_extrap=1)=0;

    virtual int int_backward(float lon,         //location of test point
		float lat,
		time_class t0,          //time of test point
		float critdist,
		time_class &tbackward,
		long &nchan,            //total number of channels
		float *val[4],          //values at eight interpolation neighbours
		float weight[4],       //their weights
		int no_extrap=1)=0;

    virtual float colloc_forward(float lon,             //location of test point
                float lat,
                time_class t0,          //time of test point
                float critdist,
                time_class &tforward,
                long &scanline,          //scan number
                long &pixnum,            //pixel number
                char *&filename,         //(pointer to) filename
                long &nchan,            //total number of channels
                float *&val)=0;            //(pointer to) values of nearest neighbour

    virtual float colloc_backward(float lon,            //location of test point
                float lat,
                time_class t0,          //time of test point
                float critdist,
                time_class &tbackward,
                long &scanline,          //scan number
                long &pixnum,            //pixel number
                char *&filename,         //(pointer to) filename
                long &nchan,            //total number of channels
                float *&val)=0;            //(pointer to) values of nearest neighbour


};

