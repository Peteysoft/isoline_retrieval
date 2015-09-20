#ifndef READ_AMSU_H
#define READ_AMSU_H 1

#include "time_class.h"

#define AMSU_A_NP 30
#define AMSU_B_NP 90

#define AMSU_A_NCHAN 15
#define AMSU_B_NCHAN 5


struct amsu_1c_rec {
  float * lon;		//longitudes in scan line
  float * lat;		//latitudes in scan line
  float **bt;		//brightness temperatures
  time_class date;	//date of scan
};

struct amsu_1c_data {
  amsu_1c_rec *data;
  long nscan;		//number of scan lines
  long np;		//points per scan
  long nchan;		//number of channels
  float missing;	//values for missing data
};

void swap_endian(unsigned long *array, long n);

long read_amsu_1c_head(char *filename, time_class & t1, time_class & t2, 
		long &np, long &nchan, int swap_end=0);

amsu_1c_data * read_amsu_1c(char *filename, int swap_end=0);

amsu_1c_data * interpolate_amsu_1c_scan(amsu_1c_data *data1, long npnew, long offset=0);
amsu_1c_data * interpolate_amsu_1c_scan2(amsu_1c_data *data1, amsu_1c_data *data2);
amsu_1c_data * cp_amsu_1c_date(amsu_1c_data *data1);
void interpolate_amsu_1c_date(amsu_1c_data *data1, amsu_1c_data *data2);

void delete_amsu_1c_data(amsu_1c_data *data);

//*** note: the 'swap_end' parameter in the following routines exists for
//          for compatibility reasons only and currently doesn't do anything
//
long rawread_amsu_1c_head(char *filename, time_class & t1, time_class & t2, 
		long &np, long &nchan, int swap_end=0);

amsu_1c_data * rawread_amsu_1c(char *filename, int swap_end=0);
int rawwrite_amsu_1c(char *filename, amsu_1c_data *data, int swap_end=0);

float ***read_amsu_lcc(char *filename, long np, long nchan);

amsu_1c_data * amsu_lc(amsu_1c_data *data, float ***c);

void delete_amsu_lcc(float ***c, long nc);

#endif

