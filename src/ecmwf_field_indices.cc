#include <stdio.h>
#include <stdlib.h>

#include "time_class.h"

#define NLON 240
#define NLAT 121

int main(int argc, char **argv) {
  unsigned long tl;		//time in longword format
  int yy, mm, dd, hh;
  time_class t0;

  char tstring[30];

  if (argc != 2) {
    printf("Syntax: ecmwf_field_indices date\n");
    printf("	where:\n");
    printf("date1:	date in format YYYYMMDDHH\n");
    exit(1);
  }

  sscanf(argv[1], "%d", &tl);
  hh=tl % 100;
  tl=tl/100;
  dd=tl % 100;
  tl=tl/100;
  mm=tl % 100;
  yy=tl/100;

  t0.init(yy, mm, dd, hh, 0, 0);
  t0.write_string(tstring);

  for (long ilat=0; ilat<NLAT; ilat++) {
    for (long ilon=0; ilon<NLON; ilon++) {
      printf("%15s%5d%5d\n", tstring, ilon, ilat);
    }
  }

}

