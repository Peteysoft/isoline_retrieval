#include <stdio.h>

#include "time_class.h"

#define NLEVELS 60
#define ECMWF_PROF_RECLEN sizeof(ecmwf_profile)

struct ecmwf_level {
  float t;		//temperature
  float q;		//water-vapour
  float cl;		//cloud liquid water
  float ci;		//cloud ice water content
  float cc;		//cloud cover
};

struct ecmwf_profile {
  //date:
  time_class date;
  //location:
  float lon;
  float lat;

  //surface variables:
  float p0;
  float t0;

  //2 m variables:
  float p2m;			//2m pressure level
  float t2m;			//temperature
  float q2m;			//water-vapour
  float u2m;			//zonal wind
  float v2m;			//meridional wind

  ecmwf_level profile[NLEVELS];

};

long write_ecmwf_profile(FILE *fs, ecmwf_profile *prof);

long read_ecmwf_profile(FILE *fs, ecmwf_profile *prof, long n);

float calc_pressure_level(ecmwf_profile *prof, long sigma);

