#include <stdio.h>

#include "read_amsu.h"

#define COMMENT_LENGTH 500

//***note: the 'swap_end' parameter exists for compatibility reasons only
//	   and currently doesn't do anything

long rawread_amsu_1c_head(char *filename, time_class & t1, time_class & t2,
		long &np, long &nchan, int swap_end) {
  FILE *fs;
  char comment[COMMENT_LENGTH];
  long nscan;

  //open the file:
  fs=fopen(filename, "r");
  if (fs == NULL) {
    fprintf(stderr, "File, %s, could not be opened for reading\n", filename);
    np=0; nchan=0;
    return -1;
  }

  //read the comment header:
  fread(comment, sizeof(char), COMMENT_LENGTH, fs);
  fprintf(stderr, "%s", comment);

  //read the attributes:
  fread(&t1, sizeof(time_class), 1, fs);
  fread(&t2, sizeof(time_class), 1, fs);
  fread(&nscan, sizeof(nscan), 1, fs);
  fread(&np, sizeof(np), 1, fs);
  fread(&nchan, sizeof(nchan), 1, fs);
  //fread(&missing, sizeof(missing), 1, fs);
 
  fclose(fs);

  return nscan;

}

int rawwrite_amsu_1c(char * filename, amsu_1c_data *data, int swap_end) {
  FILE *fs;
  long nwrit;
  char comment[COMMENT_LENGTH];
  char tstring1[30], tstring2[30];

  //open the file:
  fs=fopen(filename, "w");
  if (fs == NULL) {
    fprintf(stderr, "File, %s, could not be opened for writing\n", filename);
    return 1;
  }

  //write the headers:
 
  //comment header:
  data->data[0].date.write_string(tstring1);
  data->data[data->nscan-1].date.write_string(tstring2);
  sprintf(comment, "AMSU or AMSU-derived data in raw binary format\nDate range:%30s-%30s\n%10d scan lines; %4d points per line; %3d channels\n", 
		tstring1, tstring2, data->nscan, data->np, data->nchan);
  fwrite(comment, sizeof(char), COMMENT_LENGTH, fs);

  //attributes:
  fwrite(&data->data[0].date, sizeof(time_class), 1, fs);
  fwrite(&data->data[data->nscan-1].date, sizeof(time_class), 1, fs);
  fwrite(&data->nscan, sizeof(data->nscan), 1, fs);
  fwrite(&data->np, sizeof(data->np), 1, fs);
  fwrite(&data->nchan, sizeof(data->nchan), 1, fs);
  fwrite(&data->missing, sizeof(data->missing), 1, fs);

  //write out the data scan-line by scan-line:
  for (long i=0; i<data->nscan; i++) {
    fwrite(&data->data[i].date, sizeof(time_class), 1, fs);
    fwrite(data->data[i].lon, sizeof(float), data->np, fs);
    fwrite(data->data[i].lat, sizeof(float), data->np, fs);
    for (long j=0; j<data->np; j++) {
      fwrite(data->data[i].bt[j], sizeof(float), data->nchan, fs);
    }
  }

  fclose(fs);

  return 0;

}
      
amsu_1c_data * rawread_amsu_1c(char * filename, int swap_end) {
  FILE *fs;
  time_class t1, t2;
  char comment[COMMENT_LENGTH];
  amsu_1c_data *data;

  //open the file:
  fs=fopen(filename, "r");
  if (fs == NULL) {
    fprintf(stderr, "File, %s, could not be opened for reading\n", filename);
    return NULL;
  }
	  
  //read the comment header:
  fread(comment, sizeof(char), COMMENT_LENGTH, fs);
  fprintf(stderr, "%s", comment);

  data=new amsu_1c_data;

  //read the attributes:
  fread(&t1, sizeof(time_class), 1, fs);
  fread(&t2, sizeof(time_class), 1, fs);
  fread(&data->nscan, sizeof(data->nscan), 1, fs);
  fread(&data->np, sizeof(data->np), 1, fs);
  fread(&data->nchan, sizeof(data->nchan), 1, fs);
  fread(&data->missing, sizeof(data->missing), 1, fs);

  data->data=new amsu_1c_rec[data->nscan];

  //read the data scan-line by scan-line:
  for (long i=0; i<data->nscan; i++) {
    //allocate space:
    data->data[i].lon=new float[data->np];
    data->data[i].lat=new float[data->np];
    data->data[i].bt=new float *[data->np];

    fread(&data->data[i].date, sizeof(time_class), 1, fs);
    fread(data->data[i].lon, sizeof(float), data->np, fs);
    fread(data->data[i].lat, sizeof(float), data->np, fs);

    for (long j=0; j<data->np; j++) {
      data->data[i].bt[j]=new float[data->nchan];
      fread(data->data[i].bt[j], sizeof(float), data->nchan, fs);
    }
  }

  fclose(fs);

  return data;

}
