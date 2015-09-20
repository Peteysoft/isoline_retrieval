#include <stdio.h>
#include <stdlib.h>

#include "ecmwf_sampling2.h"

int main(int argc, char ** argv) {
  ecmwf_profile prof;
  long * slev;
  float ** vmr;		//interpolated vmr
  long nlev;
  long n;
  long mmax;

  char filename[200];
  FILE * fs;

  float mconst;		//conversion coefficient

  if (argc < 3) {
    printf("\nExtracts water-vapour volume mixing ratio at sigma\n");
    printf("levels from compressed binary files\n\n");
    printf("usage:  get_sigma_levels2 [-s] profile_file fbase level1 [level2 [level3 ...]]\n");
    printf("\nwhere:\n");
    printf("-s           = extracts specific humidity (mass-mixing-ratio)\n");
    printf("profile_file = compressed binary file containing ECMWF profiles\n");
    printf("fbase        = base name of output files\n");
    printf("level1, 2... = sigma levels to extract (as index)\n\n");
    printf("Data is output in binary form to files with the naming convention:\n\n");
    printf("<fbase><level#>s.dat\n\n");
    exit(1);
  }

  if (argv[1][0]=='-') {
    if (argv[1][1]=='s') {
      mconst=1./1.60627;
    } else {
      fprintf(stderr, "Warning: get_sigma_levels; %s unrecognized option\n", argv[1]);
    }
    argv++;
    argc--;
  } else {
    mconst=1;
  }

  nlev=argc-3;
  slev = new long[nlev];

  for (long i=0; i<nlev; i++) sscanf(argv[i+3], "%d", slev+i);

  fs = fopen(argv[1], "r");
  
  //calculate the number of records:
  fseek(fs, 0, SEEK_END);
  n=ftell(fs)/ECMWF_PROF_RECLEN;
  fseek(fs, 0, SEEK_SET);

  vmr=new float * [n];

  for (long i=0; i < n; i++) {
    read_ecmwf_profile(fs, &prof, 1);
    vmr[i]=new float [nlev];

    //interpolate in pressure:
    for (long j=0; j<nlev; j++) {
      vmr[i][j]=prof.profile[slev[j]].q;
      printf("%g ", vmr[i][j]);
    }
    printf("\n");

  }
  
  for (long i=0; i<nlev; i++) {
    sprintf(filename, "%s%2.2ds.dat", argv[2], slev[i]);
    fs=fopen(filename, "w");
    for (long j=0; j<n; j++) {
      fwrite(vmr[j]+i, 1, sizeof(float), fs);
    }
    fclose(fs);
  }

}

