#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "time_class.h"

#include "ecmwf_sampling2.h"

#define NLON 240
#define NLAT 121

#define RTTOV_NLEV 43
#define ECMWF_NLEV 60
#define NCOL 10
#define MAX_LL 200

int main(int argc, char ** argv) {
  long batchsize;
  long nbatch;
  long ibatch;

  FILE *fs1;

  long nchan;
  long instr;
  long satid;
  long isurf;

  long nvar;

  float zangle;		//zenith angle

  long ntot;		//total number of profiles

  char operand[100];	//argument of option command

  long k;

  long offset;		//batch to start at
  int head_flag;	//1 for header, 0 for none

  //for a while these variables were sitting inside the loop
  //for use on a shared-memory, parallel machine:
  ecmwf_profile *prof;
  float *bt, *bt_cld;

  char * proffile;
  char * resultsfile;

  nbatch=-1;
  offset=0;
  head_flag=1;

  FILE *infs;
  FILE *outfs;

  //parse the options:
  if (argc > 1) while(argv[1][0] == '-') {
    for (k=0; argv[1][k] != '=' && argv[1][k] != 0; k++);
    if (argv[1][k] == '=') {
      strcpy(operand, argv[1]+k+1);
      argv[1][k]=0;
      if (strcmp(argv[1], "-offset") == 0) {
        sscanf(operand, "%d", &offset);
      } else if (strcmp(argv[1], "-nbatch") == 0) {
        sscanf(operand, "%d", &nbatch);
      } else {
        printf("Unrecognized option: %s\n", argv[1]);
        exit(2);
      }
    } else {
      if (strcmp(argv[1], "-nh")==0) {
        head_flag=0;
      } else {
        printf("Unrecognized option: %s\n", argv[1]);
        exit(2);
      }
    }
    argv++;
    argc--;
  }

  if (argc != 8) {
    printf("usage:  rtsim_train8 [-offset=offset] [-nbatch=nbatch] [-nh]\n");
    printf("                     profile_file outfile sat_id \n");
    printf("                     instr_id landmask zenang batchsize\n");
    printf("	where:\n");
    printf("profile_file  = binary file listing profile locations\n");
    printf("outfile       = output binary vector file\n");
    printf("sat_id        = satellite id\n");
    printf("instr_id      = instrument id (3 = AMSU-A, 4 = AMSU-B\n");
    printf("landmask      = land-sea mask (0 = sea, 1 = land)\n");
    printf("zenangle      = zenith angle\n");
    printf("batchsize     = batch size\n");
    exit(1);
  }

  sscanf(argv[3], "%d", &satid);
  sscanf(argv[4], "%d", &instr);
  sscanf(argv[5], "%d", &isurf);
  sscanf(argv[6], "%g", &zangle);
  sscanf(argv[7], "%d", &batchsize);

  //open the input file and count the number of records:
  infs=fopen(argv[1], "r");
  fseek(infs, 0, SEEK_END);
  ntot=ftell(infs)/sizeof(ecmwf_profile);

  printf("%d records found in file %s\n", ntot, argv[1]);

  //move to offset batch in file:
  fseek(infs, offset*ECMWF_PROF_RECLEN*batchsize, SEEK_SET);
  printf("Moving to byte %d in input file\n", offset*ECMWF_PROF_RECLEN*batchsize);

  if (nbatch == -1) nbatch=(ntot-1)/batchsize-offset+1;

  printf("Preparing for %d RTTOV batches\n", nbatch);

  if (instr == 3) {
    printf("Calculating radiances for AMSU-A\n");
    nchan=15;
  } else if (instr == 4) {
    printf("Calculating radiances for AMSU-B\n");
    nchan=5;
  }

  nvar=2*nchan;

  //allocate space:
  bt=new float[nchan*batchsize];
  bt_cld=new float[nchan*batchsize];
  prof=new ecmwf_profile[batchsize];

  proffile=new char [strlen(argv[2])+14];
  resultsfile=new char [strlen(argv[2])+21];

  outfs=fopen(argv[2], "w");
  if (head_flag) fwrite(&nvar, sizeof(long), 1, outfs);

  //if (isurf ==0) se=0.95; else se=0.6;

  for (ibatch=0; ibatch<nbatch; ibatch++) {

    FILE *fs;

    long ich;
    long k, m;
    long nprof;		//number in current batch
    
    char line[MAX_LL];
    double val1, val2;

    sprintf(proffile, "%s.prof%4.4d.tmp", argv[2], ibatch);
    sprintf(resultsfile, "%s.results%4.4d.tmp", argv[2], ibatch);

    fs=fopen(proffile, "w");
    printf("Batch # %d\n", ibatch);

    nprof=read_ecmwf_profile(infs, prof, batchsize);
    printf("Number of profiles: %d\n", nprof);

    for (long iprof=0; iprof<nprof; iprof++) {
      //write this all out to the profile file:
      //write surface and 2m variables:
      fprintf(fs, "%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e\n",
                prof[iprof].lon,
                prof[iprof].lat,
                1.*isurf,
                prof[iprof].t0,
                prof[iprof].p0,
  		prof[iprof].t2m,
                prof[iprof].q2m,
                prof[iprof].u2m,
                prof[iprof].v2m);

      //write temperature:
      for (long j=0; j<ECMWF_NLEV/NCOL; j++) {
        for (k=0; k<NCOL; k++) fprintf(fs, "%16.6e", prof[iprof].profile[j*NCOL+k].t);
        fprintf(fs, "\n");
      }
      //write specific humidity:
      for (long j=0; j<ECMWF_NLEV/NCOL; j++) {
        for (k=0; k<NCOL; k++) fprintf(fs, "%16.6e", prof[iprof].profile[j*NCOL+k].q);
        fprintf(fs, "\n");
      }
      //write cloud-cover:
      for (long j=0; j<ECMWF_NLEV/NCOL; j++) {
        for (k=0; k<NCOL; k++) {
          m=j*NCOL+k;
          for (k=0; k<NCOL; k++) fprintf(fs, "%16.6e", prof[iprof].profile[m].cc);
        }
        fprintf(fs, "\n");
      }
      //write cloud liquid water content:
      for (long j=0; j<ECMWF_NLEV/NCOL; j++) {
        for (k=0; k<NCOL; k++) fprintf(fs, "%16.6e", prof[iprof].profile[j*NCOL+k].cl);
         fprintf(fs, "\n");
      }
      //write cloud ice water content:
      for (long j=0; j<ECMWF_NLEV/NCOL; j++) {
        for (k=0; k<NCOL; k++) fprintf(fs, "%16.6e", prof[iprof].profile[j*NCOL+k].ci);
        fprintf(fs, "\n");
      }
      //write precipitation, water and ice, which we take as zero
      //since there is no data...
      for (long j=0; j<ECMWF_NLEV/NCOL; j++) {
        for (k=0; k<NCOL; k++) fprintf(fs, "%16.6e", 0.0);
        fprintf(fs, "\n");
      }
      for (long j=0; j<ECMWF_NLEV/NCOL; j++) {
        for (k=0; k<NCOL; k++) fprintf(fs, "%16.6e", 0.0);
        fprintf(fs, "\n");
      }

    }
    fclose(fs);

    //if (nprof == 0) break;

    //run the "fast" radiative transfer simulator:
    fs=popen("/freax/storage/home/pmills/rtsim_bremen2/rttovscatt_test.out", "w");
    fprintf(fs, "%d\n", 2);			//Brightness temperatures
    printf("%d\n", 2);			//Brightness temperatures
    fprintf(fs, "%d\n", satid);			//satellite id (16)
    printf("%d\n", satid);			//satellite id (16)
    fprintf(fs, "%d\n", instr);			//instrument id
    printf("%d\n", instr);			//instrument id
    fprintf(fs, "%f\n", zangle);		//zenith angle
    printf("%f\n", zangle);		//zenith angle
    fprintf(fs, "%d\n", 1);		//ice crystal type (= aggregates)
    printf("%d\n", 1);		//ice crystal type (= aggregates)
    fprintf(fs, "%d\n", 1);	//ice effective size scheme (=Wyser)
    printf("%d\n", 1);		//ice effective size scheme (=Wyser)
    fprintf(fs, "%s\n", proffile);
    printf("%s\n", proffile);
    fprintf(fs, "%s\n", resultsfile);
    printf("%s\n", resultsfile);
    printf("%d\n", pclose(fs));
    //exit(1);
    //continue;

    fs=fopen(resultsfile, "r");
    for (long j=0; j<6; j++) {
      fgets(line, MAX_LL, fs);
      //printf("%s", line);
    }
    for (long i=0; i<nprof; i++) {
      for (long j=0; j<nchan; j++) {
        fgets(line, MAX_LL, fs);
        //printf("%s", line);
        k=j+i*nchan;
        sscanf(line, "%d %le %le", &ich, &val1, &val2);
        bt_cld[k]=val1;
        bt[k]=val2;
        //printf("%d %f %f\n", ich, val1, val2);
      }
      //fseek(outfs, ((ibatch*batchsize+i)*nchan*2+1)*4, SEEK_SET);
      fwrite(bt_cld+i*nchan, sizeof(float), nchan, outfs);
      fwrite(bt+i*nchan, sizeof(float), nchan, outfs);

      for (long j=0; j<nchan; j++) {
        k=j+i*nchan;
        printf("%g ", bt_cld[k]);
      }
      printf("     ");
      for (long j=0; j<nchan; j++) {
        k=j+i*nchan;
        printf("%g ", bt[k]);
      }
      printf("\n");

    }
    fclose(fs);

    //remove the temporary files:
    sprintf(line, "rm -f %s", proffile);
    system(line);
    sprintf(line, "rm -f %s", resultsfile);
    system(line);
  }
  printf("Currently at byte %d in input file\n", ftell(infs));

  delete [] bt;
  delete [] bt_cld;
  delete [] prof;

  delete [] proffile;
  delete [] resultsfile;

  fclose(infs);
  fclose(outfs);

}

