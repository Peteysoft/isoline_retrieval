#include <stdio.h>
#include <string.h>

#define MAXLL 1000

int main(int argc, char **argv) {
  FILE *fs1;
  FILE *fs2;

  char line1[MAXLL];
  char line2[MAXLL];

  long i;

  fs1=fopen(argv[1], "r");
  fs2=fopen(argv[2], "r");

  i=1;

  while (feof(fs1)==0 && feof(fs2)==0) {
    fgets(line1, MAXLL, fs1);
    fgets(line2, MAXLL, fs2);

    if (strcmp(line1, line2) != 0) {
      printf("Files differ on line %d: %s %s\n", i, line1, line2);
      fclose(fs1);
      fclose(fs2);
      return 1;
    }
    i++;
  }

  return 0;

}
  
