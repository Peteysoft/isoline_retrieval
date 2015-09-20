#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include "read_amsu.h"

int main(int argc, char ** argv) {
  amsu_1c_data * data;
  char tstring[30];

  char c;
  int swap_flag;

  amsu_1c_data * (*amsu_reader) (char *, int);

  amsu_reader=&read_amsu_1c;
  swap_flag=0;

  while ((c = getopt (argc, argv, "sR")) != -1) {
    switch (c) {
        case ('R'):
            amsu_reader=&rawread_amsu_1c;
            break;
        case ('s'):
            swap_flag=1;
            break;
        case ('?'):
            fprintf(stderr, "Unknown option: %c --ignored\n", optopt);
            break;
        default:
            fprintf(stderr, "Error parsing command line\n");
            exit(2);
    }
  }

  argc-=optind;
  argv+=optind;

  if (argc != 1) {
    printf("Usage:  print_amsu_1c_file filename\n");
    return 1;
  }

  data=(*amsu_reader) (argv[0], swap_flag);

  if (data == NULL) {
    printf("Read error, exiting... \n");
    return -1;
  }

  for (long i=0; i<data->nscan; i++) {
    for (long j=0; j<data->np; j++) {
      data->data[i].date.write_string(tstring);
      printf("%25s (%f, %f):", tstring, data->data[i].lon[j], data->data[i].lat[j]);
      for (long k=0; k<data->nchan; k++) {
        printf(" %g", data->data[i].bt[j][k]);
      }
      printf("\n");
    }
  }

  delete_amsu_1c_data(data);

}

