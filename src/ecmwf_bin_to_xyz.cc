#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#define NLON 240
#define NLAT 121

int main(int argc, char ** argv) {
  FILE *fs;
  long n;
  float * dataf;
  long * datal;
  char * tp;
  long k;

  float lon, lat;

  if (argc != 2) {
    printf("purpose: converts a binary file containing scalar data\n");
    printf("	representing a single ecmwf field into a set of coordinate\n");
    printf("	(lon-lat)/ ordinate triplets (x, y, z)\n");
    printf("	Sends results to standard out\n");
    printf("\n");
    printf("usage: ecmwf_bin_to_xyz file\n\n");
    exit(1);
  }

  fs=fopen(argv[1], "r");

  fseek(fs, 0, SEEK_END);
  n=ftell(fs)/4;
  fseek(fs, 0, SEEK_SET);

  assert(n == NLON*NLAT);

  if (argc == 3) {
    tp=argv[2];
  } else {
    tp=new char[2];
    tp[0]='f';
    tp[1]='\0';
  }

  printf("%d\n", n);

  k=0;
  if (tp[0] == 'l') {
    datal=new long[n];
    fread(datal, n, 4, fs);
    for (long j=0; j<NLAT; j++) {
      lat=j*180./(NLAT-1)-90.;
      for (long i=0; i<NLON; i++) {
        lon=i*360./NLON;
        printf("%f %f %d\n", lon, lat, datal[k]);
        k++;
      }
    }
  } else {
    dataf=new float[n];
    fread(dataf, n, 4, fs);
    for (long j=0; j<NLAT; j++) {
      lat=j*180./(NLAT-1)-90.;
      for (long i=0; i<NLON; i++) {
        lon=i*360./NLON;
        printf("%f %f %f\n", lon, lat, dataf[k]);
        k++;
      }
    }
  }

}  


