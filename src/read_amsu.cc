#include <stdio.h>
#include <math.h>

#include "read_amsu.h"

#define HEADLEN 512

#define MISSING -9999.99

void swap_endian(unsigned long *array, long n) {
  unsigned long byte1, byte2, byte3;
  for (long i=0; i<n; i++) {
    byte1=array[i] & 255UL;
    array[i] = array[i] >> 8;
    byte2=array[i] & 255UL;
    array[i] = array[i] >> 8;
    byte3=array[i] & 255UL;
    array[i] = array[i] >> 8;
    array[i]+=256UL*(byte3+256UL*(byte2+256UL*byte1));
  }
}

void delete_amsu_1c_data(amsu_1c_data *data) {
  long nscan=data->nscan;
  long np=data->np;

  for (long i=0; i<nscan; i++) {
    for (long j=0; j<np; j++) {
      delete [] data->data[i].bt[j];
    }
    delete [] data->data[i].lon;
    delete [] data->data[i].lat;
    delete [] data->data[i].bt;
  }

  delete [] data->data;
  delete data;

}

long read_amsu_1c_head(char *filename, time_class & t1, time_class & t2, 
		long &np, long &nchan, int swap_end) {

  FILE *fs;  

  long headrec[HEADLEN];

  long nscan;

  //char datestr1[30], datestr2[30];

  nchan=0;
  np=0;
  t1.init(0, 0, 0, 0, 0, 0);
  t2.init(0, 0, 0, 0, 0, 0);

  fs=fopen(filename, "r");
  if (fs == NULL) {
    printf("File %s could not be opened for reading\n", filename);
    return 0;
  }

  fread(headrec, sizeof(long), HEADLEN, fs);
  if (swap_end!=0) swap_endian((unsigned long *) headrec, HEADLEN);

  if ((headrec[21] <= 0) || ((headrec[21] && (~59L)) > 0)) {
    printf("Reading level 1c data...\n");
  } else {
    printf("Routine only defined for level 1c data\n");
    fclose(fs);
    return 0;
  }

  //printf("Grid index: %d\n", headrec[7]);

  if (headrec[7] == 10) {
    printf("  ... for AMSU-A\n");
    np=30;
    nchan=15;
  } else if (headrec[7] == 11) {
    printf("  ... for AMSU-B\n");
    np=90;
    nchan=5;
  } else {
    printf("Routine only defined for AMSU-A or -B\n");
    fclose(fs);
    return 0;
  }

  printf("satid is: %d\n", headrec[6]);

  //get the header data:
  t1.init(headrec[11], 1, 1, 0, 0, 0.);
  t1.add((double) headrec[12]-1);
  t1.add(((float) headrec[13])*0.001/86400);

  t2.init(headrec[15], 1, 1, 0, 0, 0.);
  t2.add((double) headrec[16]-1);
  t2.add(((float) headrec[17])*0.001/86400);

  nscan=headrec[18];

  fclose(fs);

  return nscan;

}

amsu_1c_data * read_amsu_1c(char *filename, int swap_end) {
  amsu_1c_rec * data;
  amsu_1c_data *all;

  long nscan;
  long np;
  long nchan;

  FILE *fs;  

  long headrec[HEADLEN];

  long * record;
  long reclen;

  long nread;
  long n;

  time_class date1, date2;
  time_class date1a, date2a;
  char datestr1[30], datestr2[30];
  long sec1, sec2;

  long lloffset;
  long toffset;
  long btoffset;

  nscan=0;
  nchan=0;
  np=0;

  double ddiff;

  fs=fopen(filename, "r");
  if (fs == NULL) {
    printf("File %s could not be opened for reading\n", filename);
    return NULL;
  }

  fread(headrec, sizeof(long), HEADLEN, fs);
  if (swap_end!=0) swap_endian((unsigned long *) headrec, HEADLEN);

  if ((headrec[21] <= 0) || ((headrec[21] && (~59L)) > 0)) {
    printf("Reading level 1c data...\n");
  } else {
    printf("Routine only defined for level 1c data\n");
    fclose(fs);
    return NULL;
  }

  //printf("Grid index: %d\n", headrec[7]);

  if (headrec[7] == 10) {
    printf("  ... for AMSU-A\n");
    reclen=768;
    np=AMSU_A_NP;
    nchan=AMSU_A_NCHAN;

    toffset=3;
    lloffset=25;
    btoffset=208;
  } else if (headrec[7] == 11) {
    printf("  ... for AMSU-B\n");
    reclen=1152;
    np=AMSU_B_NP;
    nchan=AMSU_B_NCHAN;

    toffset=3;
    lloffset=14;
    btoffset=557;
  } else {
    printf("Routine only defined for AMSU-A or -B\n");
    fclose(fs);
    return NULL;
  }

  printf("satid is: %d\n", headrec[6]);

  //get the header data:
  //printf("Year: %d, DOY: %d, milli-sec: %d\n", headrec[11], headrec[12], headrec[13]);
  date1.init(headrec[11], 1, 1, 0, 0, 0.);
  date1.add((double) headrec[12]-1);
  sec1=headrec[13];

  //printf("Year: %d, DOY: %d, milli-sec: %d\n", headrec[15], headrec[16], headrec[17]);
  date2.init(headrec[15], 1, 1, 0, 0, 0.);
  date2.add((double) headrec[16]-1);
  sec2=headrec[17];

  //ddiff=date2-date1;
  //printf("ddiff: %lf\n", ddiff);

  date1a=date1;
  date1a.add(((float) sec1)*0.001/86400);
  date1a.write_string(datestr1);
  date2a=date2;
  date2a.add(((float) sec2)*0.001/86400);
  date2a.write_string(datestr2);
  
  printf("Date range: %s - %s\n", datestr1, datestr2);

  nscan=headrec[18];

  data=new amsu_1c_rec[nscan];  

  //rewind the file:
  fseek(fs, 0, SEEK_SET);
  
  record=new long [reclen];
  nread=fread(record, sizeof(long), reclen, fs);

  for (n=0; n<nscan; n++) {
    nread=fread(record, sizeof(long), reclen, fs);
    if (nread != reclen) {
      printf("Warning:  file read error in %s, truncating data...\n", filename);
      nscan=n;
      break;
    }
    if (swap_end!=0) swap_endian((unsigned long *) record, reclen);

    //read in the date:
    data[n].date=date1;

    //if (record[toffset] < sec1 && n > nscan/2) {
    //check to see if the record is for the next day:
    if (record[toffset] < sec1 && sec2-sec1 < 0.) {
      data[n].date.add(1+record[toffset]*0.001/86400);
    } else {
      data[n].date.add(record[toffset]*0.001/86400);
    }
    //printf("%d\n", record[toffset]);

    //read in the lons and lats:
    data[n].lon=new float[np];
    data[n].lat=new float[np];
    for (long i=0; i<np; i++) {
      data[n].lat[i]=record[lloffset+2*i]*0.0001;
      data[n].lon[i]=record[lloffset+2*i+1]*0.0001;
    }

    //read in the brightness temperatures:
    data[n].bt=new float * [np];
    for (long i=0; i<np; i++) {
      data[n].bt[i]=new float[nchan];
      for (long j=0; j<nchan; j++) {
        data[n].bt[i][j]=record[btoffset+j+nchan*i]*0.01;
      }
    }

  }

  //clean up:
  delete [] record;
  fclose(fs);

  all=new amsu_1c_data;
  all->data=data;
  all->nscan=nscan;
  all->np=np;
  all->nchan=nchan;

  all->missing=MISSING;

  return all;

}

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
