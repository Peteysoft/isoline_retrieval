#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "ecmwf_sampling2.h"
#include "amsu_cret_var.h"

#include "agf_util.h"

int main(int argc, char **argv) {

  float **btA;
  float **btB;
  long mA, mB, nA, nB;

  long ndim=amsu_cret_nindA+amsu_cret_nindB+3;
  float bts[ndim];

  ecmwf_profile prof;
  FILE *ps;

  FILE *fs;

  time_class t0;
  char tstr[30];
  char tstr2[30];

  if (argc != 5) {
    printf("Collects brightness temperature data from selected channels\n");
    printf("of RTTOV simulations for AMSU-A and -B satellite instruments\n");
    printf("\n");
    printf("usage: rttov_sim_collect afile bfile pfile outfile\n");
    printf("\nwhere:\n");
    printf("	afile	= binary file containing AMSU-A simulation results\n");
    printf("	bfile	= binary file containing AMSU-B simulation results\n");
    printf("    pfile   = binary file containing list of profiles\n");
    printf("	outfile	= binary file containing output\n");
    printf("		  as a set of measurement vectors\n");
    exit(1);
  }

  btA=read_vecfile(argv[1], mA, nA);
  btB=read_vecfile(argv[2], mB, nB);

  assert(mA==mB);

  ps=fopen(argv[3], "r");

  fs=fopen(argv[4], "w");
  fwrite(&ndim, sizeof(ndim), 1, fs);

  for (long i=0; i<mA; i++) {
    read_ecmwf_profile(ps, &prof, 1);
    for (long j=0; j<amsu_cret_nindA; j++) bts[j]=btA[i][amsu_cret_indA[j]];
    for (long j=0; j<amsu_cret_nindB; j++) 
	    bts[j+amsu_cret_nindA]=btB[i][amsu_cret_indB[j]];
    //longitude and latitude:
    bts[ndim-3]=prof.lon;
    bts[ndim-2]=prof.lat;
    //time of year:
    t0.init(prof.date.year(), 1, 1, 0, 0, 0);
    bts[ndim-1]=prof.date.diff(t0);

    fwrite(bts, sizeof(float), ndim, fs);

    continue;
    t0.write_string(tstr);
    prof.date.write_string(tstr2);
    printf("%s %s %f\n", tstr, tstr2, bts[ndim-1]);
  }

  fclose(fs);
  fclose(ps);

}

