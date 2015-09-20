#include <stdio.h>

#include "amsu_interp_obj.h"

int main(int argc, char **argv) {
  amsu_interp_obj *aio;

  time_class t;
  float lon, lat;

  time_class tforward;
  float *val[4];
  float w[4];
  long nchan;

  char tstring[30];

  aio=new amsu_interp_obj("amsu_int.ini");

  //scan the command line:
  sscanf(argv[1], "%f", &lon);
  sscanf(argv[2], "%f", &lat);

  //keep t fixed:
  t.read_string("2003/1/3-12");

  aio->int_forward(lon, lat, t, 15., tforward, nchan, val, w);

  printf("values:\n");
  for (long i=0; i<nchan; i++) printf("%f ", val[0][i]);
  printf("\n");
  for (long i=0; i<nchan; i++) printf("%f ", val[1][i]);
  printf("\n");
  for (long i=0; i<nchan; i++) printf("%f ", val[2][i]);
  printf("\n");
  for (long i=0; i<nchan; i++) printf("%f ", val[3][i]);
  printf("\n");

  printf("weights:\n");
  printf("%f\n", w[0]);
  printf("%f\n", w[1]);
  printf("%f\n", w[2]);
  printf("%f\n\n", w[3]);

  tforward.write_string(tstring);
  printf("%s\n", tstring);

  delete aio;

}

