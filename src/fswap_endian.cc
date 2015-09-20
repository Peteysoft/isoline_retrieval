
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
  char *infile;		//input file name
  char *outfile;	//output file name

  FILE *ifs;		//input file stream
  FILE *outfs;		//output file stream

  char data[4];		//data to exchange
  char swp;		//"swap" byte

  long nr;		//bytes read
  long nw;		//bytes written

  int exit_value=0;
  
  if (argc < 2) {
    printf("\nByte swaps files with only 32-bit word-aligned data\n\n");
    printf("usage: fswap_endian infile [outfile]\n\n");
    printf("    where:\n");
    printf("infile  = input file name\n");
    printf("outfile = output file name; if not specified, writes to input\n\n");
    exit(-1);
  }

  //"parse" command line arguments:
  infile=argv[1];

    ifs=fopen(infile, "r");
    if (ifs == NULL) {
      fprintf(stderr, "ERROR: fswap_endian--cannot open file '%s' for reading\n", infile);
      exit(1);
    }

  if (argc > 2) {
    outfile=argv[2];
    outfs=fopen(outfile, "w");
    if (outfs == NULL) {
      fprintf(stderr, "ERROR: fswap_endian--cannot open file '%s' for writing\n", outfile);
      exit(1);
    }
  } else {
    outfile=argv[1];
    outfs=fopen(outfile, "r+");
    if (outfs == NULL) {
      fprintf(stderr, "ERROR: fswap_endian--cannot open file '%s' for writing\n", outfile);
      exit(1);
    }
  }


    while (1) {
      nr=fread(data, sizeof(char), 4, ifs);
      if (nr < 4) {
        if (feof(ifs)) {
          if (nr == 0) {
            break;
          } else {
            fprintf(stderr, "WARNING: file '%s' does not contain only word-aligned values\n", infile);
            fprintf(stderr, "now exiting to system...\n");
            break;
          }
        } else {
          fprintf(stderr, "ERROR: fswap_endian--error reading from file '%s' at byte %d, \n", 
			infile, ftell(ifs));
          fprintf(stderr, "now exiting to system...\n");
          exit_value=-2;
          break;
        }
      }


      //swap the bytes:
      swp=data[3];
      data[3]=data[0];
      data[0]=swp;
      swp=data[2];
      data[2]=data[1];
      data[1]=swp;

      nw=fwrite(data, sizeof(char), 4, outfs);

      if (nr != nw) {
        fprintf(stderr, "ERROR: fswap_endian--error writing to file, %s, at byte %d\n", outfile, ftell(outfs));
        fprintf(stderr, "now exiting to system...\n");
        exit_value=-3;
        break;
      }

    }

    fclose(ifs);
    fclose(outfs);



  return exit_value;

}

