#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "agf_util.h"

#define MAXLL 200

int main(int argc, char **argv) {
  char *infile;
  char *outfile;

  FILE *fs;		//input file stream (indices of points to remove)

  float **data;		//the vector data to prune
  float **datap;
  long m, n;		//dimensions of vector data
  long ind;		//a single index to prune

  char line[MAXLL];	//line read in
  long nless;		//number of points fewer

  if (argc != 3) {
    printf("Reads a set of indices from standard input\n");
    printf("and removes the corresponding records (vectors)\n");
    printf("from a binary file containing vector data\n");
    printf("\n");
    printf("prune_vecfile infile outfile < indices.txt");
    exit(1);
  }

  infile=argv[1];
  outfile=argv[2];

  data=read_vecfile(infile, m, n);
  datap=data;

  fs=stdin;

  nless=0;
  while (feof(fs) == 0) {

    if (fgets(line, MAXLL, fs)==NULL) break;
    if (strlen(line) == 0) break;

    sscanf(line, "%d", &ind);
    if (ind < 0 || ind >= m) continue;
    datap[ind]=NULL;
    nless++;
  }

  fs=fopen(argv[2], "w");
  fwrite(&n, sizeof(n), 1, fs);

  for (long i=0; i<m; i++) if (datap[i]!=NULL) {
    fwrite(data[i], sizeof(float), n, fs);
  }

  delete [] data[0];
  delete [] data;
  fclose(fs);

}

