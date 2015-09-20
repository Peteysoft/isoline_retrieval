#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "agf_util.h"

#define MAXLL 200

int main(int argc, char **argv) {
  char *infile;
  char *outfile;

  FILE *fs;		//input file stream (indices of points to remove)

  long *data;		//the class data to prune
  char *flags;		//which to keep and which to reject
  long n;		//number of elements in data
  long ind;		//a single index to prune

  char line[MAXLL];	//line read in
  long nless;		//number of points fewer

  if (argc != 3) {
    printf("Reads a set of indices from standard input\n");
    printf("and removes the corresponding elements from a binary file\n");
    printf("containing classification data\n");
    printf("\n");
    printf("prune_clsfile infile outfile < indices.txt");
    exit(1);
  }

  infile=argv[1];
  outfile=argv[2];

  data=read_clsfile(infile, n);
  flags=new char[n];
  for (long i=0; i<n; i++) flags[i]=1;

  fs=stdin;

  nless=0;
  while (feof(fs) == 0) {

    if (fgets(line, MAXLL, fs)==NULL) break;
    if (strlen(line) == 0) break;

    sscanf(line, "%d", &ind);
    if (ind < 0 || ind >= n) continue;
    flags[ind]=0;
    nless++;
  }

  fs=fopen(argv[2], "w");

  for (long i=0; i<n; i++) if (flags[i]) {
    fwrite(data+i, sizeof(long), 1, fs);
  }

  delete [] data;
  fclose(fs);

}

