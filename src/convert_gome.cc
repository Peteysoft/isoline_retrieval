#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "getopt.h"
#include "read_gome.h"
#include "read_spec_file.cc"

#define MAXLL 100

//#define MAXSPEC 1000

int main(int argc, char ** argv) {
  char *infile, *outfile;
  FILE *spin;
  float **spec;
  int *nchan;
  int nb;
  char line[MAXLL];
  gome_data * data;
  char tstring[30];
  int get_pmd=0;
  int def_spec=0;
  int rawread=0;
  char c;

  int b1b=1, b2b=1, b3=1, b4=1;

  while ((c = getopt (argc, argv, "vRPDb:")) != -1) {
    switch (c) {
        case ('v'):
            rg_logfile=stdout;
            break;
        case ('P'):
            get_pmd=1;
            break;
	case ('D'):
            def_spec=1;
            break;
	case ('R'):
            rawread=1;
            break;
	case ('b'):
            if (strlen(optarg) < 4) {
              fprintf(stderr, "argument for -b must have at least four characters\n");
              exit(2);
            }
            if (optarg[0]!='y') b1b=0;
            if (optarg[1]!='y') b2b=0;
            if (optarg[2]!='y') b3=0;
            if (optarg[3]!='y') b4=0;
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

  if (argc != 2 && argc != 1) {
    printf("Usage:  convert_gome [-b ????] [-D] [-P] [spectra] outfile\n");
    printf("  where:\n");
    printf("spectra   file containing desired spectra\n");
    printf("outfile   output file\n");
    printf("  options:\n");
    printf("-b ????   selected bands, choose y to include, n to exclude\n");
    printf("            where 1b, 2b, 3 and 4 are the bands; default is all four\n");
    printf("-R        file is in raw, binary format (to be added...)\n");
    printf("-D        use default spectra: 1000 channels between 290 and 790 nm\n");
    printf("            in geometric series\n");
    printf("-P        include PMD data\n");
    printf("-v        verbose\n");
    printf("\n");
    return 1;
  }

  //printf("%d\n", get_pmd);

  if (argc==2) {
    spec=read_spec_file(argv[0], nchan, nb);
    outfile=argv[1];
  } else {
    outfile=argv[0];
    spec=NULL;
    nchan=NULL;
    nb=0;
    if (def_spec) {
      spec=new float *[1];
      nchan=new int [1];
      spec[0]=default_spectra(LAMBDA0, LAMBDAN, NLAMBDA);
      nchan[0]=NLAMBDA;
      nb=1;
    }
  }

  rg_select_band(b1b, b2b, b3, b4);
  data=read_gome(stdin, spec, nchan, nb, get_pmd);

  if (data == NULL) {
    printf("Read error, exiting... \n");
    return -1;
  }

  rawwrite_gome(outfile, data);

  if (nb!=0) {
    delete [] nchan;
    for (int i=0; i<nb; i++) delete [] spec[i];
    delete [] spec;
  }

  delete_gome_data(data);

}

