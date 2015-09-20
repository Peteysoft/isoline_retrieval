#include <stdio.h>
#include <stdlib.h>

//dirt simple application to convert floating point data to a pair
//of classes based on a threshold:

int main(int argc, char ** argv) {

  FILE *fs;
  long n;
  float * data;
  float thresh;
  long * cls;

  if (argc != 4) {
    printf("Usage:  float_to_class float_file thresh class_file\n");
    printf("	where:\n");
    printf("float_file	file containing floating-point data\n");
    printf("vmr_thresh	threshold value defining classes\n");
    printf("class_file	output file containing class data\n");
    exit(1);
  }

  sscanf(argv[2], "%g", &thresh);
  fs=fopen(argv[1], "r");
  fseek(fs, 0, SEEK_END);
  n=ftell(fs)/4;
  fseek(fs, 0, SEEK_SET);
  data=new float[n];
  cls=new long[n];
  fread(data, n, sizeof(float), fs);
  fclose(fs);

  for (long i=0; i<n; i++) cls[i]=(long) (data[i] > thresh);

  fs=fopen(argv[3], "w");
  fwrite(cls, n, sizeof(long), fs);
  fclose(fs);

  delete [] data;
  delete [] cls;

}
    

