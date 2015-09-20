#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

int main(int argc, char ** argv) {
  FILE *fs;
  long n1, n2, n;
  unsigned long * data1;
  unsigned long * data2;
  unsigned long * dataf;
  int * landmask;
  long nland;
  long j, k;

  if (argc != 5) {
    printf("Unites two sets of retrieval data: one for land\n");
    printf("and one for ocean using a supplied landmask\n");
    printf("\nusage: unite_landsea_bin landfile seafile landmask outfile\n");
    printf("\nwhere:\n");
    printf("landfile = binary file containing retrievals over land\n");
    printf("seafile  = binary file containing retrievals over sea\n");
    printf("landmask = ascii file containing land-mask\n");
    printf("outfile  = binary output file\n");
    exit(1);
  }
  
  fs = fopen(argv[1], "r");
  fseek(fs, 0, SEEK_END);
  n1=ftell(fs)/4;
  fseek(fs, 0, SEEK_SET);
  data1=new unsigned long [n1];
  fread(data1, sizeof(long), n1, fs);
  fclose(fs);

  fs = fopen(argv[2], "r");
  fseek(fs, 0, SEEK_END);
  n2=ftell(fs)/4;
  fseek(fs, 0, SEEK_SET);
  data1=new unsigned long [n2];
  fread(data1, sizeof(long), n2, fs);
  fclose(fs);

  n=n1+n2;
  fs = fopen(argv[3], "r");
  landmask=new int [n];
  for (long i=0; i<n; i++) {
    fscanf(fs, "%d", landmask+i);
    nland+=landmask[i];
  }
  fclose(fs);
  assert(nland == n1);

  dataf=new unsigned long[n];
  j=0;
  k=0;
  for (long i=0; i<n; i++) {
    if (landmask[i]==1) {
      dataf[i]=data1[j];
      j++;
    } else {
      dataf[i]=data2[k];
      k++;
    }
  }

  fs=fopen(argv[4], "w");
  fwrite(dataf, sizeof(long), n, fs);
  fclose(fs);

  delete [] data1;
  delete [] data2;
  delete [] dataf;

}

