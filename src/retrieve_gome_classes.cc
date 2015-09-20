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
#include "read_spec_file.cc"
#include "read_gome.h"

#define MAXLLEN 1000
#define SWAP_ENDIAN 1

#define MISSING_RET -3.

#define NCHAN_RET 1

//#define ASCII_OUT

int main(int argc, char ** argv) {
  char *initfile;
  char *infile;
  char *outfile;
  char *specfile;

  FILE *fs;
  char line[MAXLLEN];
  char filename[MAXLLEN];

  char normfile[MAXLLEN];
  long slen;

  //amsu data:
  gome_data *data;

  //classification data:
  float **border;
  float **gradient;
  nel_ta nsamp;
  float *ave;
  float *std;

  dim_ta nvar;
  long NVAR;

  //the results:
  float *vec;
  float r;

  gome_data *results;

  int missflag;		//to search for missing data

  int err;		//error code

  char tstring[30];

  char c;
  int ncum;

  //the spectrum:
  int *nchan;
  int nb;
  float **spec;

  //the default spectrum:
  float *def_spec;

  int swap_flag;
  int raw_read;
  int get_pmd;

  raw_read=0;
  swap_flag=0;
  get_pmd=0;

  //parse the command line arguments:
  while ((c = getopt (argc, argv, "sRP")) != -1) {
    switch (c) {
        case ('R'):
	    raw_read=1;
            break;
        case ('s'):
            swap_flag=1;
            break;
        case ('P'):
            get_pmd=1;
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
    printf("Performs a classification retrieval on GOME data\n");
    printf("usage:  retrieve_gome_classes borders spectra infile outfile\n\n");
    printf("arguments:\n");
    printf("  borders      = files containing discrimination borders\n");
    printf("  spectra      = file containing desired spectra\n");
    printf("  infile       = name of input file\n");
    printf("  outfile      = name of output file\n");
    printf("\n");
    printf("options:\n");
    printf("  -s           = byte-swap input files\n");
    printf("  -P           = get PMD data\n");
    printf("  -R           = input file is in raw, binary format\n");
    printf("\n");
    return 1;
  }

  initfile=argv[0];
  specfile=argv[1];
  infile=argv[2];
  outfile=argv[3];

  //read in the spectrum:
  spec=read_spec_file(specfile, nchan, nb);
  if (get_pmd) NVAR=NPMD+NAUX; else NVAR=NAUX;
  for (long j=0; j<nb; j++) NVAR+=nchan[j];

  //if we are doing a raw-read, we need the default spectra:
  if (raw_read) {
    data=rawread_gome(infile, swap_flag);
    def_spec=default_spectra_mid(LAMBDA0, LAMBDAN, NLAMBDA);
  } else {
    fs=fopen(infile, "r");
    data=read_gome(fs, spec, nchan, nb, get_pmd);
    fclose(fs);
  }

  if (data == NULL) {
    fprintf(stderr, "Unable to open input file: %s\n", infile);
    return -2;
  }

  //read the border files:
  err=agf_read_borders(initfile, border, gradient, nsamp, nvar);
  printf("Found %d border samples: %s\n", nsamp, filename);
  if (err != 0) {
    fprintf(stderr, "Error: agf_read_borders returned error code, %d\n", err);
    fprintf(stderr, "... could not read file set, %s\n", filename);
    return err;
  }
  if (nvar != NVAR) {
    fprintf(stderr, "Border samples in %s do not have the right number (%d) of dimensions: %d\n", 
		      initfile, NVAR, nvar);
    return DIMENSION_MISMATCH;
  }

  //search for a normalization file:
  strcpy(normfile, initfile);
  strcat(normfile, ".std");
  err=read_stats(normfile, ave, std, NVAR);
  if (err != 0) {
    delete [] ave;
    delete [] std;
    ave=NULL;
    std=NULL;
  } else {
    print_stats(stdout, ave, std, NVAR);
    //normalize the border samples:
    norm_vec(border, NVAR, nsamp, ave, std);
    for (long k=0; k<NVAR; k++) {
      for (long j=0; j<nsamp; j++) gradient[j][k]/=std[k];
    }
  }
    
  //begin the retrieval:
  vec=new float[NVAR];

  results=new gome_data;
  results->nchan=nvar;
  results->npix=data->npix;
  results->data=new gome_rec[data->npix];
  results->missing=data->missing;

  //perform classification retrievals for each point:
  for (long isc=0; isc<results->npix; isc++) {
    results->data[isc].counts=new float [nvar];

    //copy the longitudes and latitudes:
    results->data[isc].lon=data->data[isc].lon;
    results->data[isc].lat=data->data[isc].lat;

    //get the measurement vector:
    if (raw_read) {
      for (int j=0; j<NAUX; j++) vec[j]=data->data[isc].counts[j];
      if (get_pmd) {
        for (int j=0; j<NPMD; j++) vec[NAUX+j]=data->data[isc].counts[j+NAUX];
        ncum=NAUX+NPMD;
      } else {
        ncum=NAUX;
      }
      for (int j=0; j<nb; j++) {
        rg_interpolate_raw(def_spec, data->data[isc].counts+NAUX+NPMD,
                        NLAMBDA, spec[j], vec+ncum, nchan[j]);
        ncum+=nchan[j];
      }
    } else {
      for (long j=0; j<NVAR; j++) {
        vec[j]=data->data[isc].counts[j];
      }
    }

    results->data[isc].counts=new float[NCHAN_RET];

    missflag=0;
    for (long i=0; i<6; i++) {
      if (vec[i] == data->missing) {
        missflag=1;
        break;
      }
    }
    if (missflag) {
      results->data[isc].counts[0]=results->missing;
      continue;
    }

#ifdef ASCII_OUT
      results->data[isc].date.write_string(tstring);
      printf("%s (%10.2f, %10.2f): (", tstring, bdata->data[isc].lon, bdata->data[isc].lat);
      for (long k=0; k<NVAR; k++) printf("%8.2f, ", vec[k]);
#endif

    //normalize the vector:
    if (ave!=NULL) norm_vec(&vec, NVAR, 1, ave, std);
    //perform the classification:
    results->data[isc].counts[0]=border_classify(border, gradient, NVAR, nsamp, vec);

#ifdef ASCII_OUT
      printf(")->%8.2f\n", results->data[isc].counts[0]);
      printf("\n");
#endif

  }

  //write the results to a file:
  if (rawwrite_gome(outfile, results) != 0) exit(UNABLE_TO_OPEN_FILE_FOR_WRITING);

  //clean up:
  delete_gome_data(data);

  delete [] border[0];
  delete [] border;
  delete [] gradient[0];
  delete [] gradient;

  if (ave != NULL) delete [] ave;
  if (std != NULL) delete [] std;

  delete [] ave;
  delete [] std;

  delete [] vec;

}


