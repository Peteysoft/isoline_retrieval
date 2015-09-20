#include <stdio.h>

#include "read_amsu.h"

int main(int argc, char **argv) {
  char *datafile;
  char *cfile;
  char *outfile;

  amsu_1c_data *data;
  float *** coeff;
  amsu_1c_data *result;

  int swap_flag;
  amsu_1c_data * (*amsu_reader) (char *, int);

  char c;
  int err;

  //set defaults:
  amsu_reader=&read_amsu_1c;
  swap_flag=0;

  //parse the command line arguments:
  while ((c = getopt (argc, argv, "rs")) != -1) {
    switch (c) {
        case ('r'):
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

  if (argc != 3) {
    printf("Purpose: applies limb-correction coefficients to AMSU bts\n");
    printf("\n");
    printf("Usage: amsu_lc [-r -c] datafile cfile outfile\n");
    printf("\n");
    printf("arguments:\n");
    printf("  datafile:   Binary file containing AMSU BTs\n");
    printf("  cfile:      ASCII coefficient file\n");
    printf("  outfile:    \"Raw\" binary output file\n");
    printf("\n");
    printf("options:\n");
    printf("  -r          Data file is in \"raw\" format\n");
    printf("  -s          Byte-swap the data file\n");
    printf("\n");
    return -1;
  }
  
  datafile=argv[0];
  cfile=argv[1];
  outfile=argv[2];

  data=(* amsu_reader) (datafile, swap_flag);
  if (data == NULL) {
    fprintf(stderr, "Unable to read input file: %s\n", datafile);
    return -2;
  }

  coeff=read_amsu_lcc(cfile, data->np, data->nchan);
  if (coeff == NULL) {
    fprintf(stderr, "Unable to open input file: %s\n", cfile);
    return -3;
  }

  result=amsu_lc(data, coeff);

  err=rawwrite_amsu_1c(outfile, result);
  if (err != 0) {
    fprintf(stderr, "Unable to write output file: %s\n", outfile);
    return err;
  }

  return 0;

}
