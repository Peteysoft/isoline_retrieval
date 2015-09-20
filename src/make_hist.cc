#include <math.h>
#include <stdio.h>

#include "agf_util.h"

int main(int argc, char **argv) {

  float *data;
  long n;

  long nbin;
  float minval;
  float maxval;
  float binsize;
  long *bins;

  float ind;

  //read in the data:
  data=(float *) read_clsfile(argv[1], n);
  printf("Found %d elements in file %s\n", n, argv[1]);

  //number of bins:
  sscanf(argv[2], "%d", &nbin);

  //find the min. and max. values:
  minval=data[0];
  maxval=data[0];
  for (long i=1; i<n; i++) {
    if (data[i] < minval) {
      minval=data[i];
    } else if (data[i] > maxval) {
      maxval=data[i];
    }
  }
  if (minval <= 0) minval=1e-7;

  printf("min: %g, max: %g\n", minval, maxval);

  //calculate the binsize:
  binsize=pow(maxval/minval, 1./nbin);
  bins=new long[nbin];
  for (long i=0; i<nbin; i++) bins[i]=0;

  printf("binsize = %g\n", binsize);

  ind=log(data[0]/minval)/log(binsize);
  if (ind < 0) ind=0;
  bins[(long) ind]++;
  for (long i=1; i<n; i++) {
    //if (data[i]==data[i-1]) printf("%d: duplicate\n", i);
    ind=log(data[i]/minval)/log(binsize);
    if (ind < 0) ind=0;
    bins[(long) ind]++;
  }

  for (long i=0; i<nbin; i++) {
    printf("%g %d\n", minval*pow(binsize, i+0.5), bins[i]);
  }

  delete [] data;
  delete [] bins;

}

