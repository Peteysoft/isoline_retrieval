#include <stdio.h>

#define MAXLL 500

//#define MAXSPEC 1000

float ** read_spec_file(char *fname, int *&nchan, int &nb) {
  FILE *spin;
  float **spec;
  char line[MAXLL];

  spin=fopen(fname, "r");
  fgets(line, MAXLL, spin);
  sscanf(line, "%d", &nb);
  nchan=new int[nb];
  spec=new float * [nb];
  for (int i=0; i<nb; i++) {
    fgets(line, MAXLL, spin);
    sscanf(line, "%d", nchan+i);
    spec[i]=new float[nchan[i]];
    for (int j=0; j<nchan[i]; j++) {
      fgets(line, MAXLL, spin);
      sscanf(line, "%f", spec[i]+j);
    }
    nchan[i]--;
  }
  fclose(spin);

  return spec;

}

