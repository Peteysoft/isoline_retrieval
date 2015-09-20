#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#define NLON 240
#define NLAT 121

int main(int argc, char ** argv) {
  FILE *fs;
  long n;
  long * datal;
  long k;

  float lon, lat;

  if (argc != 2) {
    printf("purpose: converts a binary file containing class data\n");
    printf("    covering an ecmwf field into a set of coordinate\n");
    printf("    (lon-lat)/ ordinate triplets (x, y, z)\n");
    printf("    Sends results to standard out\n");
    printf("\n");
    printf("usage: ecmwf_bin_to_xyz file\n\n");
    exit(1);
  }


  fs=fopen(argv[1], "r");

  fseek(fs, 0, SEEK_END);
  n=ftell(fs)/4;
  fseek(fs, 0, SEEK_SET);

  assert(n == NLON*NLAT);

  k=0;
  datal=new long[n];
  fread(datal, n, 4, fs);

  if (argc == 3) {
    if (argv[2][0]=='t') {
      for (long j=0; j<NLON; j++) {
        lon=j*360./NLON;
        for (long i=0; i<NLAT; i++) {
          if (datal[k] == 1) {
            lat=-90.+i*180./(NLAT-1);
            printf("%f %f\n", lon, lat);
          }
          k++;
        }
      }
      exit(0);
    }
  }
  for (long j=0; j<NLAT; j++) {
    lat=90.-j*180./(NLAT-1);
    for (long i=0; i<NLON; i++) {
      if (datal[k] == 1) {
        lon=i*360./NLON;
        printf("%f %f\n", lon, lat);
      }
      k++;
    }
  }

}  


