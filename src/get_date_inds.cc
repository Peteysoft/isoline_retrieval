#include <stdio.h>
#include <string.h>

#include "time_class.h"
#include "ecmwf_sampling2.h"

#define MAX_LL 200

int main(int argc, char **argv) {
  FILE *fs;
  char line[MAX_LL];

  time_class *tsearch;
  long nsearch;

  time_class t;

  long iline;

  ecmwf_profile prof;

  if (argc < 3) {
    printf("Purpose: extracts the record indices for each of the\n");
    printf("	profiles falling on a given set of dates.\n");
    printf("\n");
    printf("usage: get_date_inds profile_file date1 [date2 [date3... ]]\n");
    exit(1);
  }

  fs=fopen(argv[1], "r");

  nsearch=argc-2;
  tsearch=new time_class[nsearch];
  for (long i=0; i<nsearch; i++) tsearch[i].read_string(argv[i+2]);

  iline=0;
  while (feof(fs)==0) {
    read_ecmwf_profile(fs, &prof, 1);

    for (long i=0; i<nsearch; i++) if (prof.date == tsearch[i]) printf("%d\n", iline);
    iline++;
  }

  delete [] tsearch;
  fclose(fs);

}

