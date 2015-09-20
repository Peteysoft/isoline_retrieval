#include <stdio.h>
#include <string.h>

#define MAXLL 1000

int main(int argc, char **argv) {
  FILE *fs;

  char line[MAXLL];

  long i;

  long ind, n;

  fs=fopen(argv[1], "r");

  n=-1;
  if (argc > 2) {
    sscanf(argv[2], "%d", &ind);
    if (argc > 3) {
      sscanf(argv[3], "%d", &n);
    }
  } else {
    ind=0;
  }

  for (i=1; i<ind; i++) {
    fgets(line, MAXLL, fs);
    if (feof(fs) != 0) {
      fclose(fs);
      return 0;
    }
  }

  i=0;
  while (feof(fs)==0) {
    fgets(line, MAXLL, fs);
    printf("%s", line);
    i++;
    if (n > 0 && i>=n) break;
  }

  fclose(fs);

  return 0;

}
  
