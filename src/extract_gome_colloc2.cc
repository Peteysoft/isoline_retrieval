#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <getopt.h>

#include "time_class.h"
#include "read_gome.h"

#include "read_spec_file.cc"

#define FBASE_START 6
#define FDATE_START 16
#define FBASE_LEN 66

#define DATE_LEN 25

#define DEFAULT_D 100.
#define DEFAULT_T 2.5

int main(int argc, char **argv) {
  int get_pmd=0;
  char c;
  float dcrit=DEFAULT_D;
  float tcrit=DEFAULT_T;
  char *infile, *specfile, *outfile1, *outfile2;
  char tmpfile[MAXLL];
  char convert[MAXLL], *command;
  char *cp_command;
  FILE *infs, *outfs1, *outfs2;
  FILE *pfs;
  float **spec;
  int *nchan;
  int nb;
  char line[MAXLL], fname[MAXLL];
  char tstring[DATE_LEN];
  time_class date;
  float lon, lat, dlon, dlat;
  float o3;
  float tlon, tlat;
  char basepath[MAXLL], gomefile[MAXLL];
  gome_data *data;
  int32_t ndim=0;
  int k;
  int pixnum1, pixnum2;
  float d, dt;

  float *def_spec;
  int rawread=0;
  float *counts;
  int ncum;
  int endpath;
  int missflag;
  int fpclass;

  strcpy(basepath, "/data/GOME1/V4/");
  strcpy(convert, "gdp01_ex -bnynyyynnnn ");
  strcpy(tmpfile, "/tmp/int_data/temp");

  while ((c = getopt (argc, argv, "Pr:d:p:R")) != -1) {
    switch (c) {
        case ('r'):
              sscanf(optarg, "%f", &dcrit);
              break;
        case ('d'):
              sscanf(optarg, "%f", &tcrit);
              break;
        case ('P'):
            get_pmd=1;
            break;
        case ('R'):
            rawread=1;
            break;
      case ('p'):
              strcpy(basepath, optarg);
              break;
        case ('?'):
            fprintf(stderr, "Unknown option: %c --ignored\n", optopt);
            break;
        default:
            fprintf(stderr, "Error parsing command line\n");
            exit(2);
    }
  }

  argc=argc-optind;
  argv=argv+optind;

  if (argc != 3) {
    printf("\n");
    printf("usage: extract_gome_colloc2 [-P] [-r dist] [-d time] colloc spec outbase\n");
    printf("\narguments:\n");
    printf("  colloc       file containing list of collocations\n");
    printf("  spec         file containing desired spectra\n");
    printf("  outbase      base name of output files (.txt and .vec)\n");
    printf("\noptions:\n");
    //printf("  -r           data in raw, binary format\n");
    printf("  -P           get PMD data\n");
    printf("  -d dist      maximum distance [km] (default=%6.1f)\n", DEFAULT_D);
    printf("  -t time      maximum separation in time [days] (default=%5.1f)\n", DEFAULT_T);
    printf("\n");
    exit(1);
  }

  infile=argv[0];
  specfile=argv[1];
  
  outfile1=new char[strlen(argv[2])+5];
  outfile2=new char[strlen(argv[2])+5];

  strcpy(outfile1, argv[2]);
  strcpy(outfile2, argv[2]);

  strcat(outfile1, ".vec");
  strcat(outfile2, ".txt");

  spec=read_spec_file(specfile, nchan, nb);

  infs=fopen(infile, "r");

  k=0;

  outfs1=fopen(outfile1, "w");
  fwrite(&ndim, sizeof(ndim), 1, outfs1);
  outfs2=fopen(outfile2, "w");

  //all these strings are fixed length....
  cp_command=new char[strlen(basepath)+FBASE_LEN+strlen(tmpfile)+30];
  command=new char[strlen(convert)+50+strlen(tmpfile)+FBASE_LEN];

  while (feof(infs) != 1) {
    fgets(line, MAXLL, infs);
    if (feof(infs)==1) break;
    //printf(line);
    sscanf(line+DATE_LEN, "%f %f %f", &lon, &lat, &o3);
    //printf("%f %f %f\n", lon, lat, o3);

    fgets(fname, MAXLL, infs);
    if (feof(infs)==1) break;
    //printf(line);

    //read in the pixel data:
    fgets(line, MAXLL, infs);
    //printf(line);
    sscanf(line, "%d %d %f %f %f %f", &pixnum1, &pixnum2, &d, &dt, &tlon, &tlat);
    //printf("%d %d %f %f\n", pixnum1, pixnum2, d, dt);

    if (d > dcrit || dt > tcrit || o3 == 0.) continue;

    //extract the name of the original GOME lvl 0 file:
    for (int j=0; fname[j]!=0; j++) {
      if (fname[j] == '/') endpath=j;
    }

    if (strncmp(fname+endpath+1, "coords", FBASE_START) != 0) {
      fprintf(stderr, "Base filename does not match default: %6c vs %s\n", fname+endpath, "coords");
    }

    strcpy(gomefile, basepath);
    strncat(gomefile, fname+endpath+FBASE_START+FDATE_START+1, 4);
    strcat(gomefile, "/");
    strncat(gomefile, fname+endpath+FBASE_START+FDATE_START+5, 2);
    strcat(gomefile, "/");
    strncat(gomefile, fname+endpath+FBASE_START+FDATE_START+7, 2);
    strcat(gomefile, "/");
    strncat(gomefile, fname+endpath+FBASE_START+1, FBASE_LEN);

    sprintf(cp_command, "cp %s %s", gomefile, tmpfile);
    //printf("%s\n", cp_command);
    //system(cp_command);
    sprintf(command, "%s -p %5d %5d %s %s", convert, pixnum1, pixnum1, gomefile, tmpfile);
    printf("%s\n", command);
    system(command);

    sprintf(command, "%s.el1", tmpfile);
    pfs=fopen(command, "r");
    if (pfs==NULL) {
      fprintf(stderr, "Read error, command failed: %s\n", command);
      continue;
    }

    //get the gome data:
    data=read_gome(pfs, spec, nchan, nb, get_pmd);
    fclose(pfs);
    if (data==NULL) {
      fprintf(stderr, "Read error, read routine failed: %s\n", gomefile);
      continue;
    }

    if (data->npix == 0) {
      fprintf(stderr, "Read error, no data: %s\n", gomefile);
      delete_gome_data(data);
      continue;
    }

    if (pixnum1 != (int) data->data[0].counts[PIXNUMIND]) {
      fprintf(stderr, "Pixel #%d not found in file %s\n", pixnum1, gomefile);
      delete_gome_data(data);
      continue;
    }

    ndim=data->nchan;
    counts=data->data[0].counts;

    missflag=0;
    for (int j=0; j<ndim; j++) {
      if (counts[j] == data->missing) {
        missflag=1;
        break;
      }
      fpclass=fpclassify(counts[j]);
      if (fpclass == FP_INFINITE || fpclass == FP_NAN) {
        fprintf(stderr, "Warning: extract_gome_colloc detected %f\n", counts[j]);
        fprintf(stderr, "in file %s, pixel number %d\n", gomefile, pixnum1);
        missflag=1;
        break;
      }
    }
    if (missflag) {
      delete_gome_data(data);
      continue;
    }

    //write the vector data:
    fwrite(counts, sizeof(float), ndim, outfs1);

    //write the coordinate and ozone data:
    data->data[0].date.write_string(tstring);
    fprintf(outfs2, "%24s %10.4f %10.4f %14.7g %9.3f %8.3f %10.4f %10.4f %10.4f %10.4f\n", 
		tstring, data->data[0].lon, data->data[0].lat, o3, d, dt, lon, lat, tlon, tlat);
    //printf("%24s %10.4f %10.4f %14.7g %9.3f %8.3f: ", 
	//	tstring, data->data[0].lon, data->data[0].lat, o3, d, dt);
    //for (int j=0; j<ndim; j++) printf("%f ", counts[j]);
    //printf("\n");

    delete_gome_data(data);

    k++;
  }

  rewind(outfs1);

  fwrite(&ndim, sizeof(ndim), 1, outfs1);

  fclose(infs);
  fclose(outfs1);
  fclose(outfs2);

  delete [] cp_command;

  if (rawread) {
    delete [] counts;
    delete [] def_spec;
  }

}

