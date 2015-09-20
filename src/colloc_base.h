#include <math.h>

#include "read_gome.h"

#define NORMAL 0
#define RAW 1

#define BIN_NLON 360
#define BIN_NLAT 180

//base class (doesn't actually do anything...)
//contains those data and methods common to both flavours
class colloc_base {

  protected:
    //list of files and data:
    int32_t nfiles;		//total number of files
    int nloaded;		//number loaded
    int maxload;		//maximum number of loaded files allowed
    int offset;			//zero point of loaded data

    char **filelist;		//list of files in order of start date
    time_class *tstart;	//start date of each file
    time_class *tend;		//end date of each file
    gome_data **amsu_data;	//amsu data

    //file-swapping algorithm:
    gome_data *get_data_forward(long index);
    gome_data *get_data_backward(long index);

    //read routines:
    int32_t swap_end_flag;
    int4_t (*amsu_head_reader) (char *, time_class &, time_class &, int4_t &, int);
    gome_data * (*amsu_reader) (char *, int);

    char *basepath;
    char *fname;

  public:
    virtual float colloc_forward(float lon,             //location of test point
                float lat,
                time_class t0,          //time of test point
                float critdist,
                time_class &tforward,
                long &pixnum,            //pixel number
                char *&filename,         //(pointer to) filename
                long &nchan,            //total number of channels
                float *&val)=0;            //(pointer to) values of nearest neighbour

    virtual float colloc_backward(float lon,            //location of test point
                float lat,
                time_class t0,          //time of test point
                float critdist,
                time_class &tbackward,
                long &pixnum,            //pixel number
                char *&filename,         //(pointer to) filename
                long &nchan,            //total number of channels
                float *&val)=0;            //(pointer to) values of nearest neighbour


};

