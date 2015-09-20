#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

int main(int argc, char ** argv) {

  FILE *fs;

  long n, n1;
  long nland;
  int *landmask;

  long *landcls;
  long *seacls;

  if (argc != 5) {
    printf("Unites two sets of retrieval data: one for land\n");
    printf("and one for ocean using a supplied landmask\n");
    printf("\nusage: unite_landsea_cls landfile seafile landmask outfile\n");
    printf("\nwhere:\n");
    printf("landfile = binary file containing retrievals over land\n");
    printf("seafile  = binary file containing retrievals over sea\n");
    printf("landmask = ascii file containing land-mask\n");
    printf("outfile  = binary output file\n");
    exit(1);
  }

  fs=fopen(argv[2], "r");
  fseek(fs, 0, SEEK_END);
  n=ftell(fs)/4;
  fseek(fs, 0, SEEK_SET);
  landcls=new long[n];
  fread(landcls, sizeof(long), n, fs);
  fclose(fs);

  fs=fopen(argv[3], "r");
  fseek(fs, 0, SEEK_END);
  n1=ftell(fs)/4;
  fseek(fs, 0, SEEK_SET);
  seacls=new long[n];
  fread(seacls, sizeof(long), n, fs);
  fclose(fs);

  assert(n==n1);

  fs=fopen(argv[1], "r");
  landmask=new int[n];
  nland=0;
  for (long i=0; i<n; i++) {
    fscanf(fs, "%d", landmask+i);
    nland+=landmask[i];
  }
  fclose(fs);

  fs=fopen(argv[4], "w");

  for (long i=0; i<n; i++) {
    if (landmask[i] == 1) {
      fwrite(landcls+i, sizeof(long), 1, fs);
    } else {
      fwrite(seacls+i, sizeof(long), 1, fs);
    }
  }

  fclose(fs);

  delete [] landmask;
  delete [] landcls;
  delete [] seacls;
}


