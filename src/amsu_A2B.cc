#include <stdio.h>

#include "read_amsu.h"

int main(long argc, char **argv) {
  amsu_1c_data *amsua;
  amsu_1c_data *amsub;
  amsu_1c_data *amsua1;
  amsu_1c_data *amsua2;
  char tstring[30];

  amsua=read_amsu_1c(argv[1], 1);
  amsub=read_amsu_1c(argv[2], 1);

  amsua1=interpolate_amsu_1c_scan(amsua, 90, -1);
  amsua2=cp_amsu_1c_date(amsub);
  interpolate_amsu_1c_date(amsua1, amsua2);

  for (long i=0; i<amsua2->nscan; i++) {
    for (long j=0; j<amsua2->np; j++) {
      amsua2->data[i].date.write_string(tstring);
      printf("%25s %f %f", tstring, amsua2->data[i].lon[j], amsua2->data[i].lat[j]);
      for (long k=0; k<amsua2->nchan; k++) {
        printf(" %f", amsua2->data[i].bt[j][k]);
      }
      printf("\n");
    }
  }

  delete_amsu_1c_data(amsua);
  delete_amsu_1c_data(amsub);
  delete_amsu_1c_data(amsua1);
  delete_amsu_1c_data(amsua2);
  
}
