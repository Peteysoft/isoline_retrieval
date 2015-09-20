#include <stdio.h>

//we use datasets to contain the landmask:
#include "simple_temp.h"
#include "dependent_temp.h"
#include "composite_dataset.h"

//classification object:
#include "cb_cbc_obj.h"

//general numerical libraries:
#include "peteys_tmpl_lib.h"

//I/O routines:
#include "read_amsu.h"

#define MAXLLEN 1000
#define SWAP_ENDIAN 1

#define MISSING_RET -3.

int main(int argc, char ** argv) {

  FILE *fs;
  char line[MAXLLEN];
  char filename[MAXLLEN];

  //amsu data:
  amsu_1c_data *araw;
  amsu_1c_data *a2b;
  amsu_1c_data *adata;
  amsu_1c_data *bdata;
  //long na=3, nb=3;
  //long aind[3]={6, 7, 8};
  //long bind[3]={2, 3, 4};

  //classification objects:
  cb_cbc_obj *landcls;
  cb_cbc_obj *seacls; 
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
  float vec[6];
  float lors;		//land or sea?
  long cls1, cls2;
  float pdf1[2], pdf2[2];
  double * ang_int;	//angle interpolation 
  double aint;
  long aind;
  float frac;
  amsu_1c_data *results;

  int missflag;		//to search for missing data

  int err;		//error code

  char tstring[30];

  if (argc != 4) {
    printf("Performs an \"isoline retrieval\" on level 1c AMSU data\n");
    printf("usage:  retrieve_landsea afile bfile outfile\n\n");
    printf("	where:\n");
    printf("afile        = file containing AMSU-A radiances\n");
    printf("bfile        = file containing AMSU-B radiances\n");
    printf("outfile      = name of output file\n");
    printf("** Note: a landmask and initialization file must be supplied\n");
    return 1;
  }

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
  vangle=new float[nangle];
  landcls=new cb_cbc_obj [nangle];
  seacls=new cb_cbc_obj [nangle];
  for (long i=0; i<nangle; i++) {
    fgets(line, MAXLLEN, fs);
    sscanf(line, "%f\n", vangle+i);
    fgets(filename, MAXLLEN, fs);
    filename[strlen(filename)-1]='\0';
    err=landcls[i].read(filename);
    if (err != 0) {
      fprintf(stderr, "Error: cb_cbc_obj returned error code, %d\n", err);
      fprintf(stderr, "... could not read file set, %s\n", filename);
      return err;
    }
    fgets(filename, MAXLLEN, fs);
    filename[strlen(filename)-1]='\0';
    err=seacls[i].read(filename);
    if (err != 0) {
      fprintf(stderr, "Error: cb_cbc_obj returned error code, %d\n", err);
      fprintf(stderr, "... could not read file set, %s\n", filename);
      return err;
    }
  }

  //read in the amsu data:
  araw=read_amsu_1c(argv[1], SWAP_ENDIAN);
  bdata=read_amsu_1c(argv[2], SWAP_ENDIAN);

  //interpolate the AMSU-A data to AMSU-B resolution:
  a2b=interpolate_amsu_1c_scan(araw, 90, -1);
  adata=cp_amsu_1c_date(bdata);
  interpolate_amsu_1c_date(a2b, adata);

  delete_amsu_1c_data(araw);
  delete_amsu_1c_data(a2b);

  //begin the retrieval:
  results=cp_amsu_1c_date(bdata);
  results->np=adata->np;
  results->nchan=1;
  results->missing=MISSING_RET;

  //find the interpolation coefficients for the angles:
  ang_int=new double[results->np/2];
  for (long i=0; i<results->np/2; i++) {
    ang_int[i]=interpolate(vangle, nangle, (float) (50-i*100./results->np), -1);
  }

  //perform classification retrievals for each point:
  for (long isc=0; isc<results->nscan; isc++) {
    results->data[isc].bt=new float *[adata->np];
    results->data[isc].lon=new float [adata->np];
    results->data[isc].lat=new float [adata->np];
    for (long ip=0; ip<adata->np; ip++) {
      //copy the longitudes and latitudes:
      results->data[isc].lon[ip]=bdata->data[isc].lon[ip];
      results->data[isc].lat[ip]=bdata->data[isc].lat[ip];

      //get the measurement vector:
      vec[0]=adata->data[isc].bt[ip][6];
      vec[1]=adata->data[isc].bt[ip][7];
      vec[2]=adata->data[isc].bt[ip][8];
      vec[3]=bdata->data[isc].bt[ip][2];
      vec[4]=bdata->data[isc].bt[ip][3];
      vec[5]=bdata->data[isc].bt[ip][4];

//      printf("%8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n", vec[0], vec[1], vec[2], vec[3], vec[4], vec[5]);

      results->data[isc].bt[ip]=new float[1];

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
      printf("%s %10.2f %10.2f ", tstring, bdata->data[isc].lon[ip], bdata->data[isc].lat[ip]);

      //perform two retrievals--one for each of the nearest angles:
      if (lors > 0.5) {
        cls1=landcls[aind].classify(vec, pdf1);
        cls2=landcls[aind+1].classify(vec, pdf2);
	printf("%7.2f: %7.3f, %7.2f: %7.3f (l)\n", vangle[aind], pdf1[1]-pdf1[0], vangle[aind+1], pdf2[1]-pdf2[0]);
      } else {
        cls1=seacls[aind].classify(vec, pdf1);
        cls2=seacls[aind+1].classify(vec, pdf2);
	printf("%7.2f: %7.3f, %7.2f: %7.3f (s)\n", vangle[aind], pdf1[1]-pdf1[0], vangle[aind+1], pdf2[1]-pdf2[0]);
      }

      //interpolate between pixels in scan line:
      //--we store only the difference in conditional probabilities...
      results->data[isc].bt[ip][0]=(1.-frac)*(pdf1[1]-pdf1[0])+frac*(pdf2[1]-pdf2[0]);
    }
    printf("\n");
  }

  //write the results to a file:
  rawwrite_amsu_1c(argv[3], results);

  //clean up:
  delete_amsu_1c_data(adata);
  delete_amsu_1c_data(bdata);
  delete_amsu_1c_data(results);

  delete [] ang_int;

}


