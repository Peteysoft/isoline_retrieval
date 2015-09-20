#include <stdio.h>
#include "read_gome.h"

#define MAXLL 100

//#define MAXSPEC 1000

int main(int argc, char ** argv) {
  char tstring[MAXLL];
  gome_data * data;

  if (argc != 2) {
    printf("Usage:  print_gome_raw filename\n");
    return 1;
  }

  data=rawread_gome(argv[1]);

  if (data == NULL) {
    printf("Read error, exiting... \n");
    return -1;
  }

  for (long i=0; i<data->npix; i++) {
      data->data[i].date.write_string(tstring);
      printf("%25s %f %f", tstring, data->data[i].lon, data->data[i].lat);
      for (long k=0; k<data->nchan; k++) {
        printf(" %f", data->data[i].counts[k]);
      }
      printf("\n");
  }

  delete_gome_data(data);

}

