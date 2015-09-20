#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "peteys_tmpl_lib.h"
#include "ecmwf_sampling2.h"

#define P0 100000
#define KAPPA 0.2857

int main(int argc, char ** argv) {
  ecmwf_profile prof;
  float * lev;		//output theta levels
  float ** vmr;		//interpolated vmr
  long nlev;
  long n;
  long mmax;

  float theta[NLEVELS];	//profile theta levels
  float p;		//pressure
  long sind[NLEVELS];	//for sorting the levels
  float ts[NLEVELS];	//theta levels sorted
  float qs[NLEVELS];	//the sorted humidity values

  double intind;
  double frac;
  long lind;

  char filename[200];
  FILE * fs;

  float mconst;		//conversion coefficient

  if (argc < 3) {
    printf("\nExtracts water-vapour volume mixing ratio at sigma\n");
    printf("levels from compressed binary files\n\n");
    printf("usage:  get_theta_levels [-s] profile_file fbase level1 [level2 [level3 ...]]\n");
    printf("\nwhere:\n");
    printf("-s           = extracts specific humidity (mass-mixing-ratio)\n");
    printf("profile_file = binary file containing ECMWF profiles\n");
    printf("fbase        = base name of output files\n");
    printf("level1, 2... = theta levels to extract in Kelvin\n\n");
    printf("Data is output in binary form to files with the naming convention:\n\n");
    printf("<fbase><level#>th.dat\n\n");
    exit(1);
  }

  if (argv[1][0]=='-') {
    if (argv[1][1]=='s') {
      mconst=1./1.60627;
    } else {
      fprintf(stderr, "Warning: get_theta_levels; %s unrecognized option\n", argv[1]);
    }
    argv++;
    argc--;
  } else {
    mconst=1;
  }

  nlev=argc-3;
  lev = new float[nlev];

  for (long i=0; i<nlev; i++) sscanf(argv[i+3], "%g", lev+i);

  fs = fopen(argv[1], "r");
  
  //calculate the number of records:
  fseek(fs, 0, SEEK_END);
  n=ftell(fs)/ECMWF_PROF_RECLEN;
  fseek(fs, 0, SEEK_SET);

  vmr=new float * [n];

  for (long i=0; i < n; i++) {
    read_ecmwf_profile(fs, &prof, 1);
    vmr[i]=new float [nlev];

    for (long j=0; j < NLEVELS; j++) {
      p=calc_pressure_level(&prof, j);
      theta[j]=prof.profile[j].t*pow((P0/p), KAPPA);
    }

    //sort the profile by theta values:
    heapsort(theta, sind, NLEVELS);

    for (long j=0; j<NLEVELS; j++) {
      ts[j]=theta[sind[j]];
      qs[j]=prof.profile[sind[j]].q;
    }
      
    printf("P0= %g, th0= %g, thtop= %g\n", prof.p0, ts[0], ts[NLEVELS-1]);

    //interpolate in pressure:
    for (long j=0; j<nlev; j++) {
      intind=interpolate(ts, NLEVELS, lev[j]);
      lind=long(intind);
      if (lind >= NLEVELS-1) {
        lind=NLEVELS-2;
	fprintf(stderr, "Warning: [%d] theta level too high\n", i);
	fprintf(stderr, "         min: %g, max: %g\n", ts[0], ts[NLEVELS-1]);
      }
      if (lind < 0) {
        lind=0;
	fprintf(stderr, "Warning: [%d] theta level too low\n", i);
	fprintf(stderr, "         min: %g, max: %g\n", ts[0], ts[NLEVELS-1]);
      }
      frac=intind-double(lind);
      vmr[i][j]=qs[lind]*(1-frac)+qs[lind+1]*frac;
      printf("%g ", vmr[i][j]);
    }
    printf("\n");

  }
  
  for (long i=0; i<nlev; i++) {
    sprintf(filename, "%s%3.0fth.dat", argv[2], lev[i]);
    fs=fopen(filename, "w");
    for (long j=0; j<n; j++) {
      fwrite(vmr[j]+i, 1, sizeof(float), fs);
    }
    fclose(fs);
  }

}

