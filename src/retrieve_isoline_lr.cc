#include <stdio.h>
#include <assert.h>

//we use datasets to contain the landmask:
#include "simple_temp.h"
#include "dependent_temp.h"
#include "composite_dataset.h"

//classification routines:
#include "agf_lib.h"

//general numerical libraries:
#include "peteys_tmpl_lib.h"

//I/O routines:
#include "read_amsu.h"

//retrieval parameters:
#include "amsu_cret_var.h"

#define MAXLLEN 1000
#define SWAP_ENDIAN 1

#define MISSING_RET -3.

#define NVAR 6

#define NCHAN_RET 1

int main(int argc, char ** argv) {
  char *afile;
  char *bfile;
  char *outfile;

  FILE *fs;
  char line[MAXLLEN];
  char filename[MAXLLEN];

  char normfile[MAXLLEN];
  long slen;

  //amsu data:
  amsu_1c_data *adata;
  amsu_1c_data *bdata;
  //long na=3, nb=3;
  //long aind[3]={6, 7, 8};
  //long bind[3]={2, 3, 4};
  long NVAR=amsu_cret_nindA+amsu_cret_nindB;

  //classification data:
  float ***border_land;
  float ***gradient_land;
  long *nsamp_land;
  float **ave_land;
  float **std_land;

  float ***border_sea;
  float ***gradient_sea;
  long *nsamp_sea;
  float **ave_sea;
  float **std_sea;

  long nvar;

  float *vangle;
  long nangle;

  //the landmask:
  composite_dataset read_landmask;
  simple<float> *longrid;
  simple<float> *latgrid;
  dependent<float> *landmask;

  long loc, dum;

  interpol_index xind, yind;

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

  long toffseta, toffsetb;	//there may be differences in time grid btw. A & B

  char c;
  int swap_flag;
  amsu_1c_data * (*amsu_reader) (char *, int);

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

  if (argc != 3) {
    printf("Performs an \"isoline retrieval\" on level 1c AMSU data\n");
    printf("usage:  retrieve_landsea afile bfile outfile\n\n");
    printf("arguments:\n");
    printf("  afile        = file containing AMSU-A radiances\n");
    printf("  bfile        = file containing AMSU-B radiances\n");
    printf("  outfile      = name of output file\n");
    printf("\n");
    printf("options:\n");
    printf("  -s           = byte-swap input files\n");
    printf("\n");
    printf("** Note: a landmask and initialization file must be supplied:\n");
    printf("         \"landmask.ds\" and \"isoret.init\" respectively...\n");
    return 1;
  }

  afile=argv[0];
  bfile=argv[1];
  outfile=argv[2];

  vec1=new float[NVAR];
  vec2=new float[NVAR];

  //get the landmask:
  fs=fopen("landmask.ds", "r");
  read_landmask.read(fs);
  fclose(fs);

  loc=read_landmask.search_var("lon", dum);
  longrid=(simple<float> *) read_landmask.get_var(loc);
  loc=read_landmask.search_var("lat", dum);
  latgrid=(simple<float> *) read_landmask.get_var(loc);
  loc=read_landmask.search_var("landmask", dum);
  landmask=(dependent<float> *) read_landmask.get_var(loc);

  //read the border files:
  fs=fopen("isoret.init", "r");
  fgets(line, MAXLLEN, fs);
  sscanf(line, "%d", &nangle);

  //allocate the variables:
  vangle=new float[nangle];

  nsamp_land=new long[nangle];
  border_land=new float **[nangle];
  gradient_land=new float **[nangle];
  ave_land=new float *[nangle];
  std_land=new float *[nangle];

  nsamp_sea=new long[nangle];
  border_sea=new float **[nangle];
  gradient_sea=new float **[nangle];
  ave_sea=new float *[nangle];
  std_sea=new float *[nangle];

  for (long i=0; i<nangle; i++) {
    fgets(line, MAXLLEN, fs);
    sscanf(line, "%f\n", vangle+i);
    fgets(filename, MAXLLEN, fs);
    slen=strlen(filename);
    filename[slen-1]='\0';
    err=agf_read_borders(filename, border_land[i], gradient_land[i], nsamp_land[i], nvar);
    printf("Found %d border samples: %s\n", nsamp_land[i], filename);
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
    ave_land[i]=new float[NVAR];
    std_land[i]=new float[NVAR];
    err=read_stats(normfile, ave_land[i], std_land[i], NVAR);
    if (err != 0) {
      delete [] ave_land[i];
      delete [] std_land[i];
      ave_land[i]=NULL;
      std_land[i]=NULL;
    } else {
      print_stats(stdout, ave_land[i], std_land[i], NVAR);
      //normalize the border samples:
      norm_vec(border_land[i], NVAR, nsamp_land[i], ave_land[i], std_land[i]);
      for (long k=0; k<NVAR; k++) {
        for (long j=0; j<nsamp_land[i]; j++) gradient_land[i][j][k]/=std_land[i][k];
      }
    }
    
    fgets(filename, MAXLLEN, fs);
    slen=strlen(filename);
    filename[slen-1]='\0';
    err=agf_read_borders(filename, border_sea[i], gradient_sea[i], nsamp_sea[i], nvar);
    printf("Found %d border samples: %s\n", nsamp_sea[i], filename);
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
    strcpy(normfile, filename);
    strcat(normfile, ".std");
    ave_sea[i]=new float[NVAR];
    std_sea[i]=new float[NVAR];
    err=read_stats(normfile, ave_sea[i], std_sea[i], NVAR);
    if (err != 0) {
      delete [] ave_sea[i];
      delete [] std_sea[i];
      ave_sea[i]=NULL;
      std_sea[i]=NULL;
    } else {
      //normalize the border samples:
      print_stats(stdout, ave_sea[i], std_sea[i], NVAR);
      norm_vec(border_sea[i], NVAR, nsamp_sea[i], ave_sea[i], std_sea[i]);
      for (long k=0; k<NVAR; k++) {
        for (long j=0; j<nsamp_sea[i]; j++) gradient_sea[i][j][k]/=std_sea[i][k];
      }
    }
  }

  fclose(fs);

  //read in the amsu data:
  adata=(*amsu_reader)(afile, swap_flag);
  if (adata == NULL) {
    fprintf(stderr, "Unable to open input file: %s\n", afile);
    return -2;
  }
  if (adata->nchan != AMSU_A_NCHAN) {
    fprintf(stderr, "Incorrect number of channels in AMSU-A file, %s: %d\n", afile, adata->nchan);
    return -3;
  }
  if (adata->np != AMSU_A_NP) {
    fprintf(stderr, "Incorrect number of scan lines in AMSU-A file, %s: %d\n", bfile, adata->np);
    return -4;
  }
  
  bdata=(*amsu_reader)(bfile, swap_flag);
  if (bdata == NULL) {
    fprintf(stderr, "Unable to open input file: %s\n", bfile);
    return -2;
  }
  if (bdata->nchan != AMSU_B_NCHAN) {
    fprintf(stderr, "Incorrect number of channels in AMSU-B file, %s: %d\n", bfile, bdata->nchan);
    return -3;
  }
  if (bdata->np != AMSU_B_NP) {
    fprintf(stderr, "Incorrect number of scan lines in AMSU-B file, %s: %d\n", bfile, bdata->np);
    return -4;
  }

  //begin the retrieval:
  results=cp_amsu_1c_date(adata);
  results->np=adata->np;
  results->nchan=NCHAN_RET;
  results->missing=MISSING_RET;

  //find the interpolation coefficients for the angles:
  ang_int=new double[results->np/2];
  for (long i=0; i<results->np/2; i++) {
    ang_int[i]=interpolate(vangle, nangle, (float) (48.875-i*97.75/(results->np-1.)), -1);
  }

  //align the start dates of the two files, if necessary:
  if (adata->data[0].date < bdata->data[0].date) {
    for (toffseta=0; adata->data[toffseta].date<bdata->data[0].date; toffseta++);
  } else {
    toffseta=0;
  }
  if (adata->data[toffseta].date > bdata->data[0].date) {
    for (toffsetb=0; bdata->data[toffsetb].date<adata->data[toffseta].date; toffsetb++);
  } else {
    toffsetb=0;
  }

  printf("toffseta=%d; toffsetb=%d\n", toffseta, toffsetb);

  //perform classification retrievals for each point:
  for (long isc=toffseta; isc<results->nscan; isc++) {
    results->data[isc].bt=new float *[results->np];
    results->data[isc].lon=new float [results->np];
    results->data[isc].lat=new float [results->np];
    for (long ip=0; ip<adata->np; ip++) {
      //copy the longitudes and latitudes:
      results->data[isc].lon[ip]=adata->data[isc].lon[ip];
      results->data[isc].lat[ip]=adata->data[isc].lat[ip];

      //get the measurement vector:
      for (long j=0; j<amsu_cret_nindA; j++) {
        vec1[j]=adata->data[isc].bt[ip][amsu_cret_indA[j]];
      }
      for (long j=0; j<amsu_cret_nindB; j++) {
        vec1[j+amsu_cret_nindA]=bdata->data[isc*3+toffsetb].bt[ip*3+1][amsu_cret_indB[j]];
      }

      //we need two copies:
      for (long j=0; j<NVAR; j++) vec2[j]=vec1[j];

      //printf("%8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n", vec1[0], vec1[1], vec1[2], vec1[3], vec1[4], vec1[5]);

      results->data[isc].bt[ip]=new float[NCHAN_RET];

      missflag=0;
      for (long i=0; i<6; i++) {
        if (vec1[i] <= bdata->missing) {
          missflag=1;
          break;
        }
      }
      if (missflag) {
        results->data[isc].bt[ip][0]=results->missing;
	continue;
      }

      //get the landmask:
      xind=longrid->interp(bdata->data[isc].lon[ip]);
      yind=latgrid->interp(bdata->data[isc].lat[ip]);
      landmask->interpol(lors, xind, yind);

      //determine the interpolation coefficients:
      if (ip < results->np/2) aint=ang_int[ip];
		else aint=ang_int[results->np-ip-1];
      aind=(long) aint;
      if (aind >= nangle-1) {
        aind=nangle-2;
      }
      frac=aint-(double) aind;

      results->data[isc].date.write_string(tstring);
      //printf("%s %10.2f %10.2f ", tstring, bdata->data[isc].lon[ip], bdata->data[isc].lat[ip]);

      //perform two retrievals--one for each of the nearest angles:
      if (lors > 0.5) {
        //normalize the vector:
        if (ave_land[aind]!=NULL) norm_vec(&vec1, NVAR, 1, ave_land[aind], std_land[aind]);
	//perform the classification:
        r1=border_classify(border_land[aind], gradient_land[aind], NVAR, nsamp_land[aind], vec1);

	//same thing for the next nearest pointing angle:
        if (ave_land[aind+1]!=NULL) norm_vec(&vec2, NVAR, 1, ave_land[aind+1], std_land[aind+1]);
        r2=border_classify(border_land[aind+1], gradient_land[aind+1], NVAR, nsamp_land[aind+1], vec2);
	//printf("%7.2f: %7.3f, %7.2f: %7.3f (l)\n", vangle[aind], r1, vangle[aind+1], r2);
      } else {
        //normalize the vector:
        if (ave_sea[aind]!=NULL) norm_vec(&vec1, NVAR, 1, ave_sea[aind], std_sea[aind]);
	//perform the classification:
        r1=border_classify(border_sea[aind], gradient_sea[aind], NVAR, nsamp_sea[aind], vec1);

	//same thing for the next nearest pointing angle:
        if (ave_sea[aind+1]!=NULL) norm_vec(&vec2, NVAR, 1, ave_sea[aind+1], std_sea[aind+1]);
        r2=border_classify(border_sea[aind+1], gradient_sea[aind+1], NVAR, nsamp_sea[aind+1], vec2);
	//printf("%7.2f: %7.3f, %7.2f: %7.3f (s)\n", vangle[aind], r1, vangle[aind+1], r2);
      }

      //interpolate between pixels in scan line:
      //--we store only the difference in conditional probabilities...
      results->data[isc].bt[ip][0]=(1.-frac)*r1+frac*r2;
    }
    //printf("\n");
  }

  //write the results to a file:
  results->data+=toffseta;
  results->nscan-=toffseta;
  rawwrite_amsu_1c(outfile, results);
  results->data-=toffseta;

  //clean up:
  delete_amsu_1c_data(adata);
  delete_amsu_1c_data(bdata);
  delete_amsu_1c_data(results);

  for (long i=0; i<nangle; i++) {
    delete [] border_land[i][0];
    delete [] border_land[i];
    delete [] gradient_land[i][0];
    delete [] gradient_land[i];

    if (ave_land[i] != NULL) delete [] ave_land[i];
    if (std_land[i] != NULL) delete [] std_land[i];

    delete [] border_sea[i][0];
    delete [] border_sea[i];
    delete [] gradient_sea[i][0];
    delete [] gradient_sea[i];

    if (ave_sea[i] != NULL) delete [] ave_sea[i];
    if (std_sea[i] != NULL) delete [] std_sea[i];
  }

  delete [] nsamp_land;
  delete [] border_land;
  delete [] gradient_land;

  delete [] ave_land;
  delete [] std_land;

  delete [] nsamp_sea;
  delete [] border_sea;
  delete [] gradient_sea;

  delete [] ave_sea;
  delete [] std_sea;

  delete [] vangle;
  delete [] ang_int;

  delete [] vec1;
  delete [] vec2;

}


