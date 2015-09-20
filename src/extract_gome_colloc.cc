#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <getopt.h>

#include "time_class.h"
#include "read_gome.h"

#include "read_spec_file.cc"

#define FBASE_START 6
#define FBASE_LEN 22

#define DATE_LEN 25

#define DEFAULT_D 100.
#define DEFAULT_T 2.5

int main(int argc, char **argv) {
  int get_pmd=0;
  char c;
  float dcrit=DEFAULT_D;
  float tcrit=DEFAULT_T;
  char *infile, *specfile, *outfile1, *outfile2;
  FILE *infs, *outfs1, *outfs2;
  float **spec;
  int *nchan;
  int nb;
  char line[MAXLL];
  char tstring[DATE_LEN];
  time_class date;
  float lon, lat, dlon, dlat;
  float o3;
  char basepath[MAXLL], readfile[MAXLL];
  gome_data *data;
  long ndim=0;
  int k;
  int pixnum1, pixnum2;
  float d, dt;

  float *def_spec;
  int rawread=0;
  float *counts;
  int ncum;
  int endpath;

  while ((c = getopt (argc, argv, "pd:t:b:r")) != -1) {
    switch (c) {
        case ('d'):
              sscanf(optarg, "%f", &dcrit);
              break;
        case ('t'):
              sscanf(optarg, "%f", &tcrit);
              break;
        case ('p'):
            get_pmd=1;
            break;
        case ('r'):
            rawread=1;
            break;
      case ('b'):
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
    printf("usage: gome_colloc [-p] [-r dist] [-t time] colloc spec outbase\n");
    printf("\narguments:\n");
    printf("  colloc       file containing list of collocations\n");
    printf("  spec         file containing desired spectra\n");
    printf("  outbase      base name of output files (.txt and .vec)\n");
    printf("\noptions:\n");
    printf("  -r           data in raw, binary format\n");
    printf("  -p           get PMD data\n");
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

  //if we are doing a raw-read, we need the default spectra:
  if (rawread) {
    def_spec=default_spectra_mid(LAMBDA0, LAMBDAN, NLAMBDA);
    if (get_pmd) ndim=NPMD+NAUX; else ndim=NAUX;
    for (long j=0; j<nb; j++) ndim+=nchan[j];
    counts=new float[ndim];
  }

  while (feof(infs) != 1) {
    fgets(line, MAXLL, infs);
    if (feof(infs)==1) break;
    //printf(line);
    sscanf(line+DATE_LEN, "%f %f %f", &lon, &lat, &o3);
    //printf("%f %f %f\n", lon, lat, o3);

    fgets(line, MAXLL, infs);
    if (feof(infs)==1) break;
    //printf(line);

    //extract the name of the original GOME lvl 0 file:
    strcpy(readfile, basepath);
    //extract the path:
    for (int j=0; line[j]!=0; j++) {
      if (line[j] == '/') endpath=j;
    }
    strncat(readfile, line, endpath+1);
    strncat(readfile, line+FBASE_START+endpath+1, FBASE_LEN);

    //read in the pixel data:
    fgets(line, MAXLL, infs);
    //printf(line);
    sscanf(line, "%d %d %f %f", &pixnum1, &pixnum2, &d, &dt);
    //printf("%d %d %f %f\n", pixnum1, pixnum2, d, dt);

    if (d > dcrit || dt > tcrit || o3 == 0.) continue;

    if (rawread) {
      strcat(readfile, ".dat");
      data=rawread_gome(readfile);
      if (data==NULL) continue;
    } else {
      strcat(readfile, ".el1");
      //get the gome data:
      data=read_gome(readfile, spec, nchan, nb, get_pmd);
    }

    if (pixnum1 != (int) data->data[pixnum2].counts[PIXNUMIND]) {
      int fndpix=0;
      for (int j=0; j<data->npix; j++) {
        if (data->data[j].counts[PIXNUMIND] == pixnum1) {
          pixnum2=j;
          fndpix=1;
          break;
        }
      }
      if (fndpix != 1) {
        fprintf(stderr, "Pixel #%d not found in file %s\n", pixnum1, readfile);
        continue;
      }
    }

    if (rawread) {
      for (int j=0; j<NAUX; j++) counts[j]=data->data[pixnum2].counts[j];
      if (get_pmd) {
        for (int j=0; j<NPMD; j++) counts[NAUX+j]=data->data[pixnum2].counts[j+NAUX];
        ncum=NAUX+NPMD;
      } else {
        ncum=NAUX;
      }
      for (int j=0; j<nb; j++) {
        rg_interpolate_raw(def_spec, data->data[pixnum2].counts+NAUX+NPMD, 
			NLAMBDA, spec[j], counts+ncum, nchan[j]);
        ncum+=nchan[j];
      }
    } else {
      ndim=data->nchan;
      counts=data->data[pixnum2].counts;
    }

    //write the vector data:
    fwrite(counts, sizeof(float), ndim, outfs1);

    //write the coordinate and ozone data:
    data->data[pixnum2].date.write_string(tstring);
    fprintf(outfs2, "%24s %10.4f %10.4f %14.7g %9.3f %8.3f\n", 
		tstring, data->data[pixnum2].lon, data->data[pixnum2].lat, o3, d, dt);

    delete_gome_data(data);

    k++;
  }

  rewind(outfs1);

  fwrite(&ndim, sizeof(ndim), 1, outfs1);

  fclose(infs);
  fclose(outfs1);
  fclose(outfs2);

  if (rawread) {
    delete [] counts;
    delete [] def_spec;
  }

}

