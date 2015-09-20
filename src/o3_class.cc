#include <stdlib.h>
#include <string.h>

#include "agf_io.h"

#define STRMAXLEN 500

#define DATE_LEN 26

int main (int argc, char **argv) {
  char line[STRMAXLEN];

  char *infile1, *infile2;
  char *outfile;

  FILE *fs;
  cls_ta *cls;
  nel_ta n;

  float **dum;
  dim_ta nvar;

  float lon, lat, o3;
  float o3thresh;

  if (argc != 3) {
    printf("Usage: o3class basename vmrthresh\n");
    printf("where:\n");
    printf("  inbase      base name of input and output files\n");
    printf("  vmrthresh   threshold for ozone vmr\n");
    exit(-1);
  }

  sscanf(argv[2], "%f", &o3thresh);

  //lets do this the hard way:
  infile1=new char[strlen(argv[1])+5];
  strcpy(infile1, argv[1]);
  strcat(infile1, ".vec");

  infile2=new char[strlen(argv[1])+5];
  strcpy(infile2, argv[1]);
  strcat(infile2, ".txt");

  outfile=new char[strlen(argv[1])+5];
  strcpy(outfile, argv[1]);
  strcat(outfile, ".cls");

  //read the vector data to get the number of samples:
  dum=read_vecfile(infile1, n, nvar);

  fs=fopen(infile2, "r");
  cls=new cls_ta [n];
  
  for (nel_ta i=0; i<n; i++) {
    fgets(line, STRMAXLEN, fs);
    sscanf(line+DATE_LEN, "%f %f %f", &lon, &lat, &o3);
    if (o3 > o3thresh) cls[i]=1; else cls[i]=0;
  }

  fclose(fs);

  fs=fopen(outfile, "w");
  fwrite(cls, sizeof(cls_ta), n, fs);
  fclose(fs);

  delete [] cls;

}

