#include <stdio.h>
#include <stdlib.h>

int main(int argc, char ** argv) {
  FILE *fs;
  int32_t n;
  float * dataf;
  int32_t * datal;
  float val;
  char * tp;

  if (argc < 2) {
    printf("Prints the contents of a binary file containing either\n");
    printf("class or scalar floating point data to a standard out\n");
    printf("\n");
    printf("print_bin_file file [l|f]\n");
    printf("\noptions:	l for class data (longword)\n");
    printf("		f for floating point (default)\n");
    exit(1);
  }

  fs=fopen(argv[1], "r");

  fseek(fs, 0, SEEK_END);
  n=ftell(fs)/4;
  fseek(fs, 0, SEEK_SET);

  if (argc == 3) {
    tp=argv[2];
  } else {
    tp=new char[2];
    tp[0]='f';
    tp[1]='\0';
  }

  //printf("%d\n", n);

  if (tp[0] == 'l') {
    datal=new int32_t[n];
    fread(datal, n, 4, fs);
    for (int32_t i=0; i<n; i++) {
      val=(float) datal[i];
      printf("%g\n", val);
    }
  } else {
    dataf=new float[n];
    fread(dataf, n, 4, fs);
    for (int32_t i=0; i<n; i++) printf("%g\n", dataf[i]);
  }

}  


