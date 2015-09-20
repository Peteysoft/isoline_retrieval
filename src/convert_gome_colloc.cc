#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <getopt.h>

#include "time_class.h"
#include "read_gome.h"

#include "read_spec_file.cc"

#define FBASE_START 11
#define FBASE_LEN 22

#define DATE_LEN 25

#define DEFAULT_D 100.
#define DEFAULT_T 2.5

//this program should really be written in a scripting language like
//bash or perl...

int main(int argc, char **argv) {
  int get_pmd=1;
  char *infile, *specfile;
  char *path1, *path2;
  char outname1[MAXLL], outname2[MAXLL];
  FILE *infs;
  FILE *fs;
  float **spec;
  int nchan[1]={NLAMBDA};
  int nb=1;
  char line[MAXLL];
  char basepath[MAXLL], gomefile[MAXLL];
  char convert[MAXLL], command[MAXLL];
  char tempfile[MAXLL], tempfile2[MAXLL];
  gome_data *data;
  int dummy=0;
  int nooverwrite=0;
  char c;

  strcpy(basepath, "/misc/raider/GOME_LVL10/REPROC/");
  strcpy(convert, "/misc/raider/GOME_LVL10/gdp01_ex_lx -bnynyyynnnn ");

  while ((c = getopt (argc, argv, "un")) != -1) {
    switch (c) {
        case ('u'):
              dummy=1;
              break;
        case ('n'):
              nooverwrite=1;
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

  if (argc != 3 && argc != 2) {
    printf("\n");
    printf("usage: convert_gome_colloc infile [path1] path2\n");
    printf("\narguments:\n");
    printf("  infile       file containing list of collocations\n");
    printf("  path1        path for ascii output\n");
    printf("  path2        path for binary output\n");
    printf("\noptions:\n");
    printf("  -u           dummy process: print commands only\n");
    printf("  -n           existing files are not overwritten\n");
    printf("\n");
    exit(1);
  }

  if (argc == 3) {
    infile=argv[0];
    path1=argv[1];
    path2=argv[2];
  } else {
    infile=argv[0];
    path2=argv[1];
  }
  
  infs=fopen(infile, "r");

  //generate default spectra:
  spec=new float *[1];
  spec[0]=default_spectra(LAMBDA0, LAMBDAN, NLAMBDA);

  while (feof(infs) != 1) {
    fgets(line, MAXLL, infs);
    if (feof(infs)==1) break;

    fgets(line, MAXLL, infs);
    if (feof(infs)==1) break;
    //printf(line);

    //extract the name of the original GOME lvl 0 file:
    strcpy(gomefile, basepath);
    strncat(gomefile, line+FBASE_START, 4);
    strcat(gomefile, "/");
    strncat(gomefile, line+FBASE_START+4, 2);
    strcat(gomefile, "/");
    strncat(gomefile, line+FBASE_START+6, 2);
    strcat(gomefile, "/");
    strncat(gomefile, line+FBASE_START, FBASE_LEN);

    //generate the names of the output files:
    if (argc == 3) {
      strcpy(outname1, path1);
      strncat(outname1, line+FBASE_START, FBASE_LEN);
    } else {
      strcpy(outname1, "temp");
    }

    strcpy(outname2, path2);
    //remove newline:
    line[strlen(line)-1]=0;
    strcat(outname2, line);

    //if the "no-overwrite" flag is set, we test to see if the output file exists first:
    if (nooverwrite) {
      fs=fopen(outname2, "r");
      if (fs != NULL) {
        fclose(fs);
        continue;
      }
    }

    //create the command for extracting the GOME lvl 1 data:
    strcpy(command, convert);
    strcat(command, gomefile);
    strcat(command, " ");
    strcat(command, outname1);

    printf("%s\n", command);
    if (dummy != 1) if (system(command)!=0) continue;

    //read in the pixel data:
    fgets(line, MAXLL, infs);
    //printf(line);
    //printf("%d %d %f %f\n", pixnum1, pixnum2, d, dt);

    //get the gome data:
    if (dummy!=1) {
      strcat(outname1, ".el1");
      data=read_gome(outname1, spec, nchan, nb, get_pmd);

      //write it to the binary file:
      rawwrite_gome(outname2, data);

      delete_gome_data(data);
    }

  }

  fclose(infs);

  delete [] spec[0];
  delete [] spec;

}

