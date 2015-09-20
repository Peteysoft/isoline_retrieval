#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/timeb.h>

#include "gsl/gsl_rng.h"
//#include "nr.h"
#include "time_class.h"
#include "peteys_tmpl_lib.h"

#define NLON 240
#define NLAT 121

int main(int argc, char **argv) {
  gsl_rng *rann;
  long idum;

  long n;
  long lonind;
  long latind;
  long tind;
  long *ind;

  float lat;

  long nt;
  time_class tint;
  time_class tdiff;
  time_class t0, tf;
  time_class *t;

  long ntyp;
  char tstring[30];
  long i, j;
  long newind;		//index chosen to replace it

  long *sind, *indn;

  timeb now;

  if (argc != 4) {
    printf("Syntax: random_ecmwf_profiles date1 date2 n\n");
    printf("	where:\n");
    printf("date1:	first date in format YYYY/MM/DD[-HH]\n");
    printf("date1:	second date in format YYYY/MM/DD[-HH]\n");
    printf("n:		number of samples to select between the \n");
    printf("		two dates, inclusive\n");
    exit(1);
  }

  ftime(&now);
  idum=(long) now.time + ((long) now.millitm) << 7;

  rann=gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rann, idum);

  t0.read_string(argv[1]);
  tf.read_string(argv[2]);
  sscanf(argv[3], "%d", &n);

  //create the time array:
  tint.read_string("0/0/0-6");
  tdiff=tf-t0;
  nt=(long) (tdiff/tint)+1;
  t=new time_class[nt];
  t[0]=t0;
  for (i=1; i<nt; i++) {
    t[i]=t[i-1]+tint;
  }

  ind=new long[n];

  ntyp=nt*NLON*NLAT;

  //determine longitude, latitude and time indices and compose
  //them into a single longword index:
  for (i=0; i<n; i++) {
    lonind=(long) (gsl_rng_uniform(rann)*NLON);
    lat = asin(2.*gsl_rng_uniform(rann)-1.);
    latind=(long) ((0.5-lat/M_PI)*NLAT);
    tind=(long) (gsl_rng_uniform(rann)*nt);
    ind[i]=lonind+NLON*(latind+NLAT*tind);
  }

  //sort:
  heapsort_inplace(ind, n);
/*
  sind=heapsort(ind, n);
  indn=map_vector(ind, sind, n);
  delete [] ind;
  ind=indn;
*/

  //find the duplicates and replace them:
  for (i=1; i<n; i++) {
    if (ind[i] == ind[i-1]) {
      j=0;
      while (j < n) {
        lonind=(long) (gsl_rng_uniform(rann)*NLON);
        lat = asin(2*gsl_rng_uniform(rann)-1);
        latind=(long) ((0.5-lat/M_PI)*NLAT);
        tind=(long) (gsl_rng_uniform(rann)*nt);
        newind=lonind+NLON*(latind+NLAT*tind);
        for (j=0; j<n; j++) if (newind == ind[j]) break;
      }
      ind[i-1]=newind;
    }
  }

  //re-sort:
  heapsort_inplace(ind, n);
/*
  sind=heapsort(ind, n);
  indn=map_vector(ind, sind, n);
  delete [] ind;
  ind=indn;
*/

  //decompose the indices and write them out to a file:
  for (long i=0; i<n; i++) {
    tind=ind[i]/NLON/NLAT;
    latind=ind[i]/NLON-tind*NLAT;
    lonind=ind[i] % NLON;
    t[tind].write_string(tstring);
    //printf("%8d %20s %5d %5d\n", ind[i], tstring, lonind, latind);
    printf("%15s%5d%5d\n", tstring, lonind, latind);
  }

  delete [] t;
  delete [] ind;

  gsl_rng_free(rann);

}

