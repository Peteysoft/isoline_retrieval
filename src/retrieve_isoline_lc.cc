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
  char *cfile_land;
  char *cfile_sea;
  char *afile_land;
  char *afile_sea;
  char *bfile;
  char *outfile;

  FILE *fs;
  char line[MAXLLEN];
  char filename[MAXLLEN];

  char normfile[MAXLLEN];
  long slen;

  //amsu data:
  amsu_1c_data *araw;
  amsu_1c_data *a2b;
  amsu_1c_data *adata_land;
  amsu_1c_data *adata_sea;
  amsu_1c_data *bdata;
  //long na=3, nb=3;
  //long aind[3]={6, 7, 8};
  //long bind[3]={2, 3, 4};
  long NVAR=amsu_cret_nindA+amsu_cret_nindB;

  //classification data:
  float **border_land;
  float **gradient_land;
  long nsamp_land;
  float *ave_land;
  float *std_land;

  float **border_sea;
  float **gradient_sea;
  long nsamp_sea;
  float *ave_sea;
  float *std_sea;

  long nvar;

  //the landmask:
  composite_dataset read_landmask;
  simple<float> *longrid;
  simple<float> *latgrid;
  dependent<float> *landmask;

  long loc, dum;

  interpol_index xind, yind;

  //the results:
  float *vec;
  float lors;		//land or sea?
  float r;

  amsu_1c_data *results;

  int missflag;		//to search for missing data

  int err;		//error code

  char tstring[30];

  char c;
  int swap_flag;
  amsu_1c_data * (*amsu_reader) (char *, int);

  amsu_reader=&rawread_amsu_1c;
  swap_flag=0;

  if (argc != 7) {
    printf("Performs an \"isoline retrieval\" on level 1c AMSU data\n");
    printf("usage:  retrieve_isoline_lc border_land border_sea afile_land afile_sea bfile outfile\n\n");
    printf("arguments:\n");
    printf("  border_land  = base of files containing class-borders trained for land\n");
    printf("  border_sea   = base of files containing class-borders trained for sea\n");
    printf("  afile_land   = AMSU-A radiances limb-corrected for land\n");
    printf("  afile_sea    = AMSU-A radiances limb-corrected for sea\n");
    printf("  bfile        = file containing limb-corrected AMSU-B radiances\n");
    printf("  outfile      = name of output file\n");
    printf("\n");
    printf("** Note: a landmask file must be supplied: \"landmask.ds\"\n");
    return 1;
  }

  cfile_land=argv[1];
  cfile_sea=argv[2];
  afile_land=argv[3];
  afile_sea=argv[4];
  bfile=argv[5];
  outfile=argv[6];

  vec=new float[NVAR];

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
  err=agf_read_borders(cfile_land, border_land, gradient_land, nsamp_land, nvar);
  if (err != 0) {
    fprintf(stderr, "Error: agf_read_borders returned error code, %d\n", err);
    fprintf(stderr, "... could not read file set, %s\n", cfile_land);
    return err;
  }
  if (nvar != NVAR) {
      fprintf(stderr, "Border samples in %s do not have the right number (%d) of dimensions: %d\n", 
		      cfile_land, NVAR, nvar);
      return DIMENSION_MISMATCH;
  }
  printf("Found %d border samples: %s\n", nsamp_land, cfile_land);

  //search for a normalization file:
  strcpy(normfile, cfile_land);
  strcat(normfile, ".std");
  ave_land=new float[NVAR];
  std_land=new float[NVAR];
  err=read_stats(normfile, ave_land, std_land, NVAR);
  if (err != 0) {
    delete [] ave_land;
    delete [] std_land;
    ave_land=NULL;
    std_land=NULL;
  } else {
    print_stats(stdout, ave_land, std_land, NVAR);
    //normalize the border samples:
    norm_vec(border_land, NVAR, nsamp_land, ave_land, std_land);
    for (long k=0; k<NVAR; k++) {
      for (long j=0; j<nsamp_land; j++) gradient_land[j][k]/=std_land[k];
    }
  }
    
  err=agf_read_borders(cfile_sea, border_sea, gradient_sea, nsamp_sea, nvar);
  if (err != 0) {
    fprintf(stderr, "Error: agf_read_borders returned error code, %d\n", err);
    fprintf(stderr, "... could not read file set, %s\n", cfile_sea);
    return err;
  }
  if (nvar != NVAR) {
    fprintf(stderr, "Border samples in %s do not have the right number (%d) of dimensions: %d\n", 
		      cfile_sea, NVAR, nvar);
    return DIMENSION_MISMATCH;
  }
  printf("Found %d border samples: %s\n", nsamp_sea, cfile_sea);

  strcpy(normfile, cfile_sea);
  strcat(normfile, ".std");
  ave_sea=new float[NVAR];
  std_sea=new float[NVAR];
  err=read_stats(normfile, ave_sea, std_sea, NVAR);
  if (err != 0) {
    delete [] ave_sea;
    delete [] std_sea;
    ave_sea=NULL;
    std_sea=NULL;
  } else {
    //normalize the border samples:
    print_stats(stdout, ave_sea, std_sea, NVAR);
    norm_vec(border_sea, NVAR, nsamp_sea, ave_sea, std_sea);
    for (long k=0; k<NVAR; k++) {
      for (long j=0; j<nsamp_sea; j++) gradient_sea[j][k]/=std_sea[k];
    }
  }

  //read in the amsu data:
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

  araw=(*amsu_reader)(afile_land, swap_flag);
  if (araw == NULL) {
    fprintf(stderr, "Unable to open input file: %s\n", afile_land);
    return -2;
  }
  if (araw->nchan != AMSU_A_NCHAN) {
    fprintf(stderr, "Incorrect number of channels in AMSU-A file, %s: %d\n", afile_land, araw->nchan);
    return -3;
  }
  
  if (araw->np == AMSU_A_NP) {
    //interpolate the AMSU-A data to AMSU-B resolution:
    a2b=interpolate_amsu_1c_scan(araw, 90, -1);
    adata_land=cp_amsu_1c_date(bdata);
    interpolate_amsu_1c_date(a2b, adata_land);

    delete_amsu_1c_data(araw);
    delete_amsu_1c_data(a2b);
  } else if (araw->np == AMSU_B_NP) {
    adata_land=araw;
  } else {
    fprintf(stderr, "Incorrect number of scan lines in AMSU-A file, %s: %d\n", afile_land, araw->np);
    return -4;
  }

  araw=(*amsu_reader)(afile_sea, swap_flag);
  if (araw == NULL) {
    fprintf(stderr, "Unable to open input file: %s\n", afile_sea);
    return -2;
  }
  if (araw->nchan != AMSU_A_NCHAN) {
    fprintf(stderr, "Incorrect number of channels in AMSU-A file, %s: %d\n", afile_sea, araw->nchan);
    return -3;
  }

  if (araw->np == AMSU_A_NP) {
    //interpolate the AMSU-A data to AMSU-B resolution:
    a2b=interpolate_amsu_1c_scan(araw, 90, -1);
    adata_sea=cp_amsu_1c_date(bdata);
    interpolate_amsu_1c_date(a2b, adata_sea);

    delete_amsu_1c_data(araw);
    delete_amsu_1c_data(a2b);
  } else if (araw->np == AMSU_B_NP) {
    adata_sea=araw;
  } else {
    fprintf(stderr, "Incorrect number of scan lines in AMSU-A file, %s: %d\n", afile_sea, araw->np);
    return -4;
  }

  //begin the retrieval:
  results=cp_amsu_1c_date(bdata);
  results->np=bdata->np;
  results->nchan=NCHAN_RET;
  results->missing=MISSING_RET;

  //perform classification retrievals for each point:
  for (long isc=0; isc<results->nscan; isc++) {
    results->data[isc].bt=new float *[results->np];
    results->data[isc].lon=new float [results->np];
    results->data[isc].lat=new float [results->np];
    for (long ip=0; ip<results->np; ip++) {
      //copy the longitudes and latitudes:
      results->data[isc].lon[ip]=bdata->data[isc].lon[ip];
      results->data[isc].lat[ip]=bdata->data[isc].lat[ip];

      //land or sea?:
      xind=longrid->interp(bdata->data[isc].lon[ip]);
      yind=latgrid->interp(bdata->data[isc].lat[ip]);
      landmask->interpol(lors, xind, yind);

      //get the measurement vector:
      if (lors > 0.5) {
        for (long j=0; j<amsu_cret_nindA; j++) {
          vec1[j]=adata_land->data[isc].bt[ip][amsu_cret_indA[j]];
        }
      } else {
        for (long j=0; j<amsu_cret_nindA; j++) {
          vec1[j]=adata_sea->data[isc].bt[ip][amsu_cret_indA[j]];
        }
      }
      //get the measurement vector:
      for (long j=0; j<amsu_cret_nindB; j++) {
        vec[j+amsu_cret_nindA]=bdata->data[isc].bt[ip][amsu_cret_indB[j]];
      }

      vec[3]=bdata->data[isc].bt[ip][2];
      vec[4]=bdata->data[isc].bt[ip][3];
      vec[5]=bdata->data[isc].bt[ip][4];

      //printf("%8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n", vec[0], vec[1], vec[2], vec[3], vec[4], vec[5]);

      results->data[isc].bt[ip]=new float[NCHAN_RET];

      missflag=0;
      for (long i=0; i<6; i++) {
        if (vec[i] <= bdata->missing) {
          missflag=1;
          break;
        }
      }
      if (missflag) {
        results->data[isc].bt[ip][0]=results->missing;
	continue;
      }

      //results->data[isc].date.write_string(tstring);
      //printf("%s %10.2f %10.2f ", tstring, bdata->data[isc].lon[ip], bdata->data[isc].lat[ip]);

      if (lors > 0.5) {
        //normalize the vector:
        if (ave_land!=NULL) norm_vec(&vec, NVAR, 1, ave_land, std_land);
        //perform the classification:
        r=border_classify(border_land, gradient_land, NVAR, nsamp_land, vec);

      } else {
        //normalize the vector:
        if (ave_sea!=NULL) norm_vec(&vec, NVAR, 1, ave_sea, std_sea);
	//perform the classification:
        r=border_classify(border_sea, gradient_sea, NVAR, nsamp_sea, vec);

      }

      //--store the difference in conditional probabilities...
      results->data[isc].bt[ip][0]=r;
    }
    //printf("\n");
  }

  //write the results to a file:
  rawwrite_amsu_1c(outfile, results);

  //clean up:
  delete_amsu_1c_data(adata_land);
  delete_amsu_1c_data(adata_sea);
  delete_amsu_1c_data(bdata);
  delete_amsu_1c_data(results);

  delete [] border_land[0];
  delete [] border_land;
  delete [] gradient_land[0];
  delete [] gradient_land;

  if (ave_land != NULL) delete [] ave_land;
  if (std_land != NULL) delete [] std_land;

  delete [] border_sea[0];
  delete [] border_sea;
  delete [] gradient_sea[0];
  delete [] gradient_sea;

  if (ave_sea != NULL) delete [] ave_sea;
  if (std_sea != NULL) delete [] std_sea;

  delete [] vec;

}


