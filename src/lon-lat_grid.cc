#include "stdio.h"
//#include "agf_util.h"

#define LON_MIN -180.
#define LON_MAX 180.
#define LAT_MIN -90.
#define LAT_MAX 90.

int main(int argc, char **argv) {
  long D;
  long xn;
  long yn;
  float lon, lat;
  float **coords;
  FILE *fs;

  sscanf(argv[1], "%d", &xn);
  sscanf(argv[2], "%d", &yn);

  coords=new float *[xn*yn];
  coords[0]=new float[2*xn*yn];
  for (long i=0; i<xn*yn; i++) coords[i]=coords[0]+2*i;

  for (long i=0; i<xn; i++) {
    lon=LON_MIN+(float) i*(LON_MAX-LON_MIN)/xn;
    for (long j=0; j<yn; j++) {
      lat=LAT_MIN+(float) j*(LAT_MAX-LAT_MIN)/(yn-1);
      printf("%8.2f %8.2f\n", lon, lat);
      coords[j*xn+i][0]=lon;
      coords[j*xn+i][1]=lat;
    }
  }

  fs=fopen(argv[3], "w");
  D=2;
  fwrite(&D, sizeof(long), 1, fs);
  fwrite(coords[0], sizeof(float), xn*yn*2, fs);
  fclose(fs);

  delete [] coords[0];
  delete [] coords;

}
  
