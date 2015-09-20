#include <stdio.h>

#include "time_class.h"

typedef int int4_t;

struct gome_rec {
  time_class date;
  float lon;
  float lat;
  float *counts;
};

struct gome_data {
  int4_t npix;		//number of pixels
  int4_t nchan;		//number of channels
  gome_rec * data;
  float missing;
};

extern FILE * rg_logfile;

#define LAMBDA0 290.
#define LAMBDAN 790.
#define NLAMBDA 1000

float * default_spectra(float lambda0, float lambdan, int n);
float * default_spectra_mid(float lambda0, float lambdan, int n);
void rg_interpolate_raw(float *specraw, float *countsraw, int nraw, 
		float *spec, float *counts, int nchan);

#define NPMD 16*3
#define NAUX 6
#define PIXNUMIND 4

void rg_select_band(int b1b, int b2b, int b3, int b4);
gome_data * read_gome(FILE *in, float **spec, int *nchan, int nb, int get_pmd=0);
void delete_gome_data(gome_data *data);
int4_t rawread_gome_head(char *filename, time_class & t1, time_class & t2,
                int4_t &nchan, int swap_end=0);
int4_t rawwrite_gome(char * filename, gome_data *data, int swap_end=0);
gome_data * rawread_gome(char * filename, int swap_end=0);

