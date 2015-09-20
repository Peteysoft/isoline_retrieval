//general-purpose isoline retrieval program:

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include <getopt.h>

//classification routines:
#include "agf_lib.h"

//general numerical libraries:
#include "peteys_tmpl_lib.h"

//I/O routines:
#include "read_amsu.h"

#define MAXLLEN 1000
#define SWAP_ENDIAN 1

#define MISSING_RET -3.

#define NCHAN_RET 1

//#define ASCII_OUT

int main(int argc, char ** argv) {
  char *initfile;
  char *infile;
  char *outfile;

  FILE *fs;
  char line[MAXLLEN];
  char filename[MAXLLEN];

  char normfile[MAXLLEN];
  long slen;

  //amsu data:
  amsu_1c_data *data;

  //classification data:
  float ***border;
  float ***gradient;
  nel_ta *nsamp;
  float **ave;
  float **std;

  dim_ta nvar;
  long NVAR;

  float *vangle;
  long nangle;

  //the results:
  float *vec1, *vec2;
  float lors;		//land or sea?
  float r1, r2;
  double * ang_int;	//angle interpolation 
  double aint;
  long aind;
  float frac;
  amsu_1c_data *results;

  int missflag;		//to search for missing data

  int err;		//error code

  char tstring[30];

  char c;
  int swap_flag;
  amsu_1c_data * (*amsu_reader) (char *, int);

  amsu_reader=&read_amsu_1c;
  swap_flag=0;

  //parse the command line arguments:
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

  if (argc != 3) {
    printf("Performs a general-purpose \"isoline retrieval\" for a cross-track-\n");
    printf("scanning satellite instrument.\n");
    printf("usage:  retrieve_isoline_g initfile infile outfile\n\n");
    printf("arguments:\n");
    printf("  initfile     = initialization file\n");
    printf("  infile       = file containing measurement vectors\n");
    printf("  outfile      = name of output file\n");
    printf("\n");
    printf("options:\n");
    printf("  -s           = byte-swap input files\n");
    printf("\n");
    return 1;
  }

  initfile=argv[0];
  infile=argv[1];
  outfile=argv[2];

  //read in the swath data:
  data=(*amsu_reader)(infile, swap_flag);
  if (data == NULL) {
    fprintf(stderr, "Unable to open input file: %s\n", infile);
    return -2;
  }

  //for now we just take the whole measurement vector:
  NVAR=data->nchan;

  //read the border files:
  fs=fopen(initfile, "r");
  fgets(line, MAXLLEN, fs);
  sscanf(line, "%d", &nangle);

  //allocate the variables:
  vangle=new float[nangle];

  nsamp=new nel_ta[nangle];
  border=new float **[nangle];
  gradient=new float **[nangle];
  ave=new float *[nangle];
  std=new float *[nangle];

  for (long i=0; i<nangle; i++) {
    fgets(line, MAXLLEN, fs);
    sscanf(line, "%f\n", vangle+i);
    fgets(filename, MAXLLEN, fs);
    slen=strlen(filename);
    filename[slen-1]='\0';
    err=agf_read_borders(filename, border[i], gradient[i], nsamp[i], nvar);
    printf("Found %d border samples: %s\n", nsamp[i], filename);
    if (err != 0) {
      fprintf(stderr, "Error: agf_read_borders returned error code, %d\n", err);
      fprintf(stderr, "... could not read file set, %s\n", filename);
      return err;
    }
    if (nvar != NVAR) {
      fprintf(stderr, "Border samples in %s do not have the right number (%d) of dimensions: %d\n", 
		      filename, NVAR, nvar);
      return DIMENSION_MISMATCH;
    }
    //search for a normalization file:
    strcpy(normfile, filename);
    strcat(normfile, ".std");
    ave[i]=new float[NVAR];
    std[i]=new float[NVAR];
    err=read_stats(normfile, ave[i], std[i], NVAR);
    if (err != 0) {
      delete [] ave[i];
      delete [] std[i];
      ave[i]=NULL;
      std[i]=NULL;
    } else {
      print_stats(stdout, ave[i], std[i], NVAR);
      //normalize the border samples:
      norm_vec(border[i], NVAR, nsamp[i], ave[i], std[i]);
      for (long k=0; k<NVAR; k++) {
        for (long j=0; j<nsamp[i]; j++) gradient[i][j][k]/=std[i][k];
      }
    }
    
  }

  fclose(fs);

  //begin the retrieval:
  vec1=new float[NVAR];
  vec2=new float[NVAR];

  results=cp_amsu_1c_date(data);
  results->np=data->np;
  results->nchan=NCHAN_RET;
  results->missing=MISSING_RET;

  //pre-calculate the interpolation coefficients for the angles:
  ang_int=new double[results->np/2];
  for (long i=0; i<results->np/2; i++) {
    ang_int[i]=interpolate(vangle, nangle, (float) (50-i*100./(results->np-1)), -1);
  }

  //perform classification retrievals for each point:
  for (long isc=0; isc<results->nscan; isc++) {
    results->data[isc].bt=new float *[results->np];
    results->data[isc].lon=new float [results->np];
    results->data[isc].lat=new float [results->np];
    for (long ip=0; ip<data->np; ip++) {
      //copy the longitudes and latitudes:
      results->data[isc].lon[ip]=data->data[isc].lon[ip];
      results->data[isc].lat[ip]=data->data[isc].lat[ip];

      //get the measurement vector:
      for (long j=0; j<NVAR; j++) {
        vec1[j]=data->data[isc].bt[ip][j];
      }

      //we need two copies:
      for (long j=0; j<NVAR; j++) vec2[j]=vec1[j];

      results->data[isc].bt[ip]=new float[NCHAN_RET];

      missflag=0;
      for (long i=0; i<6; i++) {
        if (vec1[i] == data->missing) {
          missflag=1;
          break;
        }
      }
      if (missflag) {
        results->data[isc].bt[ip][0]=results->missing;
	continue;
      }

      //determine the interpolation coefficients:
      if (ip < results->np/2) {
        aint=ang_int[ip];
      } else {
        aint=ang_int[results->np-ip-1];
      }
      aind=(long) aint;
      if (aind >= nangle-1) {
        aind=nangle-2;
      }
      frac=aint-(double) aind;

#ifdef ASCII_OUT
      results->data[isc].date.write_string(tstring);
      printf("%s (%10.2f, %10.2f): (", tstring, bdata->data[isc].lon[ip], bdata->data[isc].lat[ip]);
      for (long k=0; k<NVAR; k++) printf("%8.2f, ", vec1[k]);
#endif

      //perform two retrievals--one for each of the nearest angles:
      //normalize the vector:
      if (ave[aind]!=NULL) norm_vec(&vec1, NVAR, 1, ave[aind], std[aind]);
      //perform the classification:
      r1=border_classify(border[aind], gradient[aind], NVAR, nsamp[aind], vec1);

      //same thing for the next nearest pointing angle:
      if (ave[aind+1]!=NULL) norm_vec(&vec2, NVAR, 1, ave[aind+1], std[aind+1]);
      r2=border_classify(border[aind+1], gradient[aind+1], NVAR, nsamp[aind+1], vec2);
      //printf("%7.2f: %7.3f, %7.2f: %7.3f (l)\n", vangle[aind], r1, vangle[aind+1], r2);

      //interpolate between pixels in scan line:
      //--we store only the difference in conditional probabilities...
      results->data[isc].bt[ip][0]=(1.-frac)*r1+frac*r2;
#ifdef ASCII_OUT
      printf(")->%8.2f\n", results->data[isc].bt[ip][0]);
#endif

    }

    //printf("\n");
#ifdef ASCII_OUT
    printf("\n");
#endif
  }

  //write the results to a file:
  if (rawwrite_amsu_1c(outfile, results) != 0) exit(UNABLE_TO_OPEN_FILE_FOR_WRITING);

  //clean up:
  delete_amsu_1c_data(data);

  for (long i=0; i<nangle; i++) {
    delete [] border[i][0];
    delete [] border[i];
    delete [] gradient[i][0];
    delete [] gradient[i];

    if (ave[i] != NULL) delete [] ave[i];
    if (std[i] != NULL) delete [] std[i];

  }

  delete [] nsamp;
  delete [] border;
  delete [] gradient;

  delete [] ave;
  delete [] std;

  delete [] vangle;
  delete [] ang_int;

  delete [] vec1;
  delete [] vec2;

}


