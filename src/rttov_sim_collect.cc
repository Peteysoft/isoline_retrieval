#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "amsu_cret_var.h"

#include "agf_util.h"

int main(int argc, char **argv) {

  float **btA;
  float **btB;
  long mA, mB, nA, nB;

  long ndim=amsu_cret_nindA+amsu_cret_nindB;
  float bts[ndim];

  FILE *fs;

  if (argc != 4) {
    printf("Collects brightness temperature data from selected channels\n");
    printf("of RTTOV simulations for AMSU-A and -B satellite instruments\n");
    printf("\n");
    printf("usage: rttov_sim_collect afile bfile outfile\n");
    printf("\nwhere:\n");
    printf("	afile	= binary file containing AMSU-A simulation results\n");
    printf("	bfile	= binary file containing AMSU-B simulation results\n");
    printf("	outfile	= binary file containing output\n");
    printf("		  as a set of measurement vectors\n");
    exit(1);
  }

  btA=read_vecfile(argv[1], mA, nA);
  btB=read_vecfile(argv[2], mB, nB);

  assert(mA==mB);

  fs=fopen(argv[3], "w");
  fwrite(&ndim, sizeof(ndim), 1, fs);

  for (long i=0; i<mA; i++) {
    for (long j=0; j<amsu_cret_nindA; j++) bts[j]=btA[i][amsu_cret_indA[j]];
    for (long j=0; j<amsu_cret_nindB; j++) 
	    bts[j+amsu_cret_nindA]=btB[i][amsu_cret_indB[j]];

    fwrite(bts, sizeof(float), ndim, fs);
  }

  fclose(fs);

}

