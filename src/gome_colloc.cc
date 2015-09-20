#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "agf_util.h"
#include "amsu_tinterp_obj.h"

#define DEFAULT_DCRIT 50
#define DEFAULT_NFILES 20

#define DEFAULT_DT 0.05
#define DEFAULT_NT 20

#define MAXLL 400

#define MISSING -2.

#define TRECLEN 24

#define NW 4

//#define ASCII_OUT

int main(int argc, char **argv) {
  long nchan;			//number of channels

  amsu_interp_base *amsu_interp;	//interpolation engine
  int amsu_read_type;
  int swap_flag;
  
  float dcrit;			//critical distance
  float d2;			//critical distance squared
  int nfiles;			//number of files to keep loaded

  float *bt1, *bt2;		//brightness temperatures
  float err1, err2;

  char *tstring;
  char tstring1[30], tstring2[30], tstring3[30];
  time_class t0;		//creates interpolated field at this time

  time_class tforward, tbackward;
  float tdiff;
  double tw1, tw2;
  FILE *fs;
  FILE *infs;

  //filenames:
  char *initfile;
  char *infile;
  char *outfile;
  char *basepath;		//base path

  char *nfile=NULL;
  char *sfile=NULL;

  char c;			//for parsing the command options

  //diagnostics:
  float maxdt;
  float mindt;
  float maxtdiff;
  float mintdiff;
  float dtave;

  int missflag;

  double onesec=1./86400.;

  interpol_index dt;
  long nt;

  char line[MAXLL];
  char *check;
 
  float lon, lat, dlon, dlat; 

  long n, k;
  long scannum1, scannum2;
  long pnum1, pnum2;

  char *fname1, *fname2;

  int use2;		//flag--use backword search, default is forward

  //set defaults:
  dcrit=DEFAULT_DCRIT;
  nfiles=DEFAULT_NFILES;
  amsu_read_type=NORMAL;
  swap_flag=0;
  basepath=NULL;

  dt=DEFAULT_DT;
  nt=DEFAULT_NT;

  //parse the command line arguments:
  while ((c = getopt (argc, argv, "d:t:s:o:b:r:n:")) != -1) {
    switch (c) {
      case ('r'):
	      sscanf(optarg, "%f", &dcrit);
	      break;
      case ('n'):
	      sscanf(optarg, "%d", &nfiles);
              break;
      case ('b'):
	      basepath=new char[strlen(optarg)+1];
	      strcpy(basepath, optarg);
              break;
      case ('o'):
	      nfile=new char[strlen(optarg)+1];
	      strcpy(nfile, optarg);
	      break;
      case ('s'):
	      sfile=new char[strlen(optarg)+1];
	      strcpy(sfile, optarg);
	      break;
      case ('d'):
	      sscanf(optarg, "%f", &dt);
	      break;
      case ('t'):
	      sscanf(optarg, "%d", &nt);
	      break;
      case ('?'):
	      fprintf(stderr, "Unknown option: %c --ignored\n", optopt);
	      break;
      default:
	      fprintf(stderr, "Error parsing command line\n");
	      exit(2);

    }
  }    

  argc=argc-optind;
  argv=argv+optind;
	
  if (argc < 3) {
    printf("\n");
    printf("usage: gome_colloc [-r dcrit] [-n nfiles] [-b basepath]\n");
    printf("                         indexfile infile outfile\n");
    printf("\narguments:\n");
    printf("  indexfile    initialization file--indexes data files\n");
    printf("  infile       input file\n");
    printf("  outfile      output file\n");
    printf("\noptions:\n");
    printf("  -r dcrit     max. distance of test point to measurement pixel in km\n");
    printf("               (default=%d)\n", DEFAULT_DCRIT);
    printf("  -n nfiles    max. number of files to keep in RAM (default=%d)\n", DEFAULT_NFILES);
    printf("  -b basepath  base directory containing AMSU l1c files\n");
    printf("\n");
    exit(1);
  }

  initfile=argv[0];
  infile=argv[1];
  outfile=argv[2];

  //check to make sure the output file is openable:
  fs=fopen(outfile, "w");
  if (fs == NULL) {
    fprintf(stderr, "Unable to open output file: %s\n", outfile);
    exit(3);
  }

  infs=fopen(infile, "r");
  
  //initialize the interpolation object:
  if (nfile != NULL && sfile != NULL) {
    amsu_interp=new amsu_tinterp_obj(initfile, nfiles, basepath, nfile, sfile, dt, nt);
  } else {
    amsu_interp=new amsu_interp_obj(initfile, nfiles, basepath);
  }

  n=0;
  while (feof(infs) != 1) {
    check=fgets(line, MAXLL, infs);
    if (check == NULL) continue;
    t0.read_string(line);
    sscanf(line+TRECLEN, "%f %f %f %f", &lon, &lat, &dlon, &dlat);

    //perform the interpolation:
    maxdt=0;
    mindt=1000;
    maxtdiff=0;
    mintdiff=1000;
    dtave=0;

    d2=dcrit*dcrit;
    err1=amsu_interp->colloc_forward(lon+dlon, lat+dlat, t0, d2, 
		    tforward, scannum1, pnum1, fname1, nchan, bt1);


    //tforward.write_string(tstring1);
    //fprintf(stdout, "Date: %s\n", tstring1);
    
    err2=amsu_interp->colloc_backward(lon+dlon, lat+dlat, t0, d2, 
		    tbackward, scannum2, pnum2, fname2, nchan, bt2);

    //tbackward.write_string(tstring1);
    //fprintf(stdout, "Date: %s\n", tstring1);
    //
    if (err1 >= 0) tw1=tforward.diff(t0);
    if (err2 >= 0) tw2=t0.diff(tbackward);

    if (err1 < 0 && err2 < 0) continue;
    use2=0;
    if ((err1 >= 0 && err2 >= 0) && tw1 > tw2) use2=1;
    if (err1 < 0) use2=1;

    if (use2) {
      pnum1=pnum2;
      scannum1=scannum2;
      bt1=bt2;
      fname1=fname2;
      err1=err2;
      tw1=-tw2;
    }
    
    tforward.write_string(tstring1);
    t0.write_string(tstring2);
    tbackward.write_string(tstring3);

    printf("(%8.2f, %8.2f) ", lon+dlon, lat+dlat);
    printf("%s  %s ", tstring1, tstring3);
    printf("%g %g %g %2g\n", tw2, tw1, err1, err2);
    
    //printf("%s: %f+%f+%f+%f=%f\n", tstring1, w1[0], w1[1], w1[2], w1[3],
//		    w1[0]+w1[1]+w1[2]+w1[3]);
    //printf("%s: %f+%f+%f+%f=%f\n", tstring2, w2[0], w2[1], w2[2], w2[3],
//		    w2[0]+w2[1]+w2[2]+w2[3]);

    //bug searching:  if the total of the weights is much greater than 1, exit...
    //assert(w1[0]+w1[1]+w1[2]+w1[3] < 1.1);
    //assert(w2[0]+w2[1]+w2[2]+w2[3] < 1.1);

    //diagnostics:
    if (tw1 > maxdt) {
      maxdt=tw1;
    } else if (tw1 < mindt) {
      mindt=tw1;	   
    }
    dtave+=tw1;

    fprintf(fs, "%s", line);
    fprintf(fs, "%s\n", fname1);
    fprintf(fs, "%6.0f %ld %ld %g %g %f %f\n", bt1[PIXNUMIND], scannum1, pnum1, sqrt(err1), tw1, $
		amsu_interp->traj_lon, amsu_interp->traj_lat);

    n++;
  }

  fclose(fs);
  fclose(infs);

  //print the diagnostics:
  dtave/=n;
  fprintf(stderr, "diagnostic parameter          %8s   %8s   %8s\n",
	                     "min", "max", "average");
  fprintf(stderr, "time to scan line:           %10.3g %10.3g %10.3g\n",
		  mindt, maxdt, dtave/2);
  fprintf(stderr, "time interval:               %10.3g %10.3g %10.3g\n\n",
		  mintdiff, maxtdiff, dtave);

  if (basepath != NULL) delete [] basepath;

}

