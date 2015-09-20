#include "read_amsu.h"

int main(int argc, char ** argv) {
  char *infile;
  char *outfile;

  amsu_1c_data *data;
  amsu_1c_data *cropped;
  amsu_1c_data * (*amsu_reader) (char *, int);
  int swap_flag;
  long ncrop;

  char c;

  amsu_reader=&read_amsu_1c;
  swap_flag=0;


  //parse the command line arguments:
  while ((c = getopt (argc, argv, "sr")) != -1) {
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

  infile=argv[0];
  sscanf(argv[1], "%d", &ncrop);
  outfile=argv[2];

  data=(* amsu_reader) (infile, swap_flag);

  cropped=cp_amsu_1c_date(data);
  cropped->np=data->np-2*ncrop;
  cropped->nchan=data->nchan;

  for (long i=0; i<cropped->nscan; i++) {
    //cropped->data[i].lon=new float [cropped->np];
    //cropped->data[i].lat=new float [cropped->np];
    //cropped->data[i].bt=new float * [cropped->np];

    for (long j=ncrop; j<data->np-ncrop; j++) {
      //cropped->data[i].bt[j-ncrop]=new float[cropped->nchan];
      cropped->data[i].bt=data->data[i].bt+ncrop;
      //cropped->data[i].lon[j-ncrop]=data->data[i].lon[j];
      cropped->data[i].lon=data->data[i].lon+ncrop;
      //cropped->data[i].lat[j-ncrop]=data->data[i].lat[j];
      cropped->data[i].lat=data->data[i].lat+ncrop;
      //for (long k=0; k<cropped->nchan; k++) cropped->data[i].bt[j-ncrop][k]=data->data[i].bt[j][k];
    }
  }

  rawwrite_amsu_1c(outfile, cropped);

  delete_amsu_1c_data(data);
  delete [] cropped->data;
  delete cropped;
  //delete_amsu_1c_data(cropped);

  return 0;

}
