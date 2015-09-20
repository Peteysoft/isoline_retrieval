#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "agf_io.h"
#include "amsu_tinterp_obj.h"

#define DEFAULT_DCRIT 50
#define DEFAULT_NFILES 20

#define DEFAULT_DT 0.05
#define DEFAULT_NT 20

#define MISSING -2.

//#define ASCII_OUT

int main(int argc, char **argv) {
  float **coords;		//lon-lat coordinates
  dim_ta n;			//number of interpolation points
  nel_ta D;			//number of dimensions (should be 2)
  long nchan;			//number of channels

  amsu_interp_obj *amsu_interp;	//interpolation engine
  int amsu_read_type;
  int swap_flag;
  
  float dcrit;			//critical distance
  float d2;			//critical distance squared
  int nfiles;			//number of files to keep loaded

  float w1[4];			//weights
  float w2[4];			//weights
  float *bt1[4], *bt2[4];	//brightness temperatures
  int err1, err2;

  float **bt_interp;		//interpolated brightness temperatures

  char *tstring;
  char tstring1[30], tstring2[30], tstring3[30];
  time_class t0;		//creates interpolated field at this time

  time_class tforward, tbackward;
  float tdiff;
  float tw1, tw2;
  FILE *fs;

  //filenames:
  char *initfile;
  char *coordfile;
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

  //set defaults:
  dcrit=DEFAULT_DCRIT;
  nfiles=DEFAULT_NFILES;
  amsu_read_type=NORMAL;
  swap_flag=0;
  basepath=NULL;

  dt=DEFAULT_DT;
  nt=DEFAULT_NT;

  //parse the command line arguments:
  while ((c = getopt (argc, argv, "h:t:S:N:p:r:n:")) != -1) {
    switch (c) {
      case ('r'):
	      sscanf(optarg, "%f", &dcrit);
	      break;
      case ('n'):
	      sscanf(optarg, "%d", &nfiles);
              break;
      case ('p'):
	      basepath=new char[strlen(optarg)+1];
	      strcpy(basepath, optarg);
              break;
      case ('N'):
	      nfile=new char[strlen(optarg)+1];
	      strcpy(nfile, optarg);
	      break;
      case ('S'):
	      sfile=new char[strlen(optarg)+1];
	      strcpy(sfile, optarg);
	      break;
      case ('h'):
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
	
  if (argc < 4) {
    printf("\n");
    printf("usage: interpolate_amsu [-r dcrit] [-n nfiles] [-p basepath]\n");
    printf("                         indexfile date coordfile outfile\n");
    printf("\narguments:\n");
    printf("  indexfile    initialization file--indexes data files\n");
    printf("  t0           interpolation time\n");
    printf("  coordfile    binary file containing list of lon-lat coordinates\n");
    printf("  outfile      binary output file\n");
    printf("\noptions:\n");
    printf("  -r dcrit     max. distance of test point to measurement pixel in km\n");
    printf("               (default=%d)\n", DEFAULT_DCRIT);
    printf("  -n nfiles    max. number of files to keep in RAM (default=%d)\n", DEFAULT_NFILES);
    printf("  -p basepath  base directory containing AMSU l1c files\n");
    printf("\n");
    exit(1);
  }

  initfile=argv[0];
  tstring=argv[1];
  coordfile=argv[2];
  outfile=argv[3];

  //read in the lon-lat pairs:
  coords=read_vecfile(coordfile, n, D);
  assert(D==2);

  //check to make sure the output file is openable:
  fs=fopen(outfile, "w");
  if (fs == NULL) {
    fprintf(stderr, "Unable to open output file: %s\n", outfile);
    exit(3);
  }
  
  //initialize the interpolation object:
  if (nfile != NULL && sfile != NULL) {
    amsu_interp=new amsu_tinterp_obj(initfile, nfiles, basepath, nfile, sfile, dt, nt);
  } else {
    amsu_interp=new amsu_interp_obj(initfile, nfiles, basepath);
  }

  t0.read_string(tstring);

  //initialize the results array:
  bt_interp=new float *[n];

  //perform the interpolation:
  maxdt=0;
  mindt=1000;
  maxtdiff=0;
  mintdiff=1000;
  dtave=0;

  d2=dcrit*dcrit;
  for (long i=0; i<n; i++) {
    tforward=t0;
    do {
      err1=amsu_interp->int_forward(coords[i][0], coords[i][1], tforward, d2, 
		    tforward, nchan, bt1, w1);
      missflag=0;

      for (long j=0; j<nchan; j++) {
        //printf("bts(1): %f %f %f %f\n", bt1[0][j], bt1[1][j], bt1[2][j], bt1[3][j]);
        if (bt1[0][j] <= MISSING || bt1[1][j] <= MISSING || bt1[2][j] <= MISSING
		    || bt1[3][j] <= MISSING ) {
          missflag=1;

          tforward.write_string(tstring1);
          tforward.add(onesec);
          tforward.write_string(tstring2);
	  fprintf(stdout, "Missing value found, %s, starting at %s\n", tstring1, tstring2);

	  break;
        }
      }
    } while (missflag);

    //tforward.write_string(tstring1);
    //fprintf(stdout, "Date: %s\n", tstring1);
    
    if (err1==1) {
      fprintf(stderr, "Coordinate (%f, %f) not covered forward by data\n",
		      coords[i][0], coords[i][1]);
      fprintf(stderr, "Please add extra entries to initialization file\n");
      exit(4);
    } else if (err1 == -1) {
      fprintf(stderr, "Warning: doing an extrapolation for (%8.2f, %8.2f): %s\n", 
		      coords[i][0], coords[i][1], tstring1);
      fprintf(stderr, "Weights: %10.3g %10.3g %10.3g %10.3g\n", w1[0], w1[1], w1[2], w1[3]);
    }

    //calculate the interpolates:
    bt_interp[i]=new float[nchan];
    for (long j=0;j<nchan; j++) {
      bt_interp[i][j]=(bt1[0][j]*w1[0]+bt1[1][j]*w1[1]+
		      bt1[2][j]*w1[2]+bt1[3][j]*w1[3])/
		      (w1[0]+w1[1]+w1[2]+w1[3]);
    }

    tbackward=t0;
    do {
      err2=amsu_interp->int_backward(coords[i][0], coords[i][1], tbackward, d2, 
		    tbackward, nchan, bt2, w2);
      missflag=0;
      for (long j=0; j<nchan; j++) {
        //printf("bts(2): %f %f %f %f\n", bt2[0][j], bt2[1][j], bt2[2][j], bt2[3][j]);
        if (bt2[0][j] == MISSING || bt2[1][j] == MISSING || bt2[2][j] == MISSING
		    || bt2[3][j] == MISSING ) {
          missflag=1;

          tbackward.write_string(tstring1);
          tbackward.add(-onesec);
          tbackward.write_string(tstring2);
	  fprintf(stdout, "Missing value found, %s, starting at %s\n", tstring1, tstring2);

	  break;
        }
      }
    } while (missflag);

    //tbackward.write_string(tstring1);
    //fprintf(stdout, "Date: %s\n", tstring1);
    
    if (err2==1) {
      fprintf(stderr, "Coordinate (%f, %f) not covered backward by data\n",
		      coords[i][0], coords[i][1]);
      fprintf(stderr, "Please add extra entries to initialization file\n");
      exit(4);
    } else if (err2 == -1) {
      fprintf(stderr, "Warning: doing an extrapolation for (%8.2f, %8.2f): %s\n", 
		      coords[i][0], coords[i][1], tstring3);
      fprintf(stderr, "Weights: %10.3g %10.3g %10.3g %10.3g\n", w2[0], w2[1], w2[2], w2[3]);
    }
    
    tforward.write_string(tstring1);
    t0.write_string(tstring2);
    tbackward.write_string(tstring3);

    //calculate weights in time:
    tdiff=tforward.diff(tbackward);
    tw1=t0.diff(tbackward)/tdiff;
    tw2=tforward.diff(t0)/tdiff;

    //calculate the interpolates:
    for (long j=0;j<nchan; j++) {
      bt_interp[i][j]=tw1*bt_interp[i][j]+tw2*(bt2[0][j]*w2[0]+bt2[1][j]*w2[1]+
		      bt2[2][j]*w2[2]+bt2[3][j]*w2[3])/
	      		(w2[0]+w2[1]+w2[2]+w2[3]);
    }

    printf("(%8.2f, %8.2f) ", coords[i][0], coords[i][1]);
    printf("%s  %s ", tstring1, tstring3);
    printf("%7.3f %7.3f %7.3f %2d %2d\n", tw2*tdiff, tw1*tdiff, tdiff, err1, err2);
    
    //printf("%s: %f+%f+%f+%f=%f\n", tstring1, w1[0], w1[1], w1[2], w1[3],
//		    w1[0]+w1[1]+w1[2]+w1[3]);
    //printf("%s: %f+%f+%f+%f=%f\n", tstring2, w2[0], w2[1], w2[2], w2[3],
//		    w2[0]+w2[1]+w2[2]+w2[3]);

    //since we don't want to extrapolate to the test point,
    //but we're too lazy to ensure that it always
    //falls within the scan-line...
    /*
    if (w1[0]<0) w1[0]=0;
    if (w1[1]<0) w1[1]=0;
    if (w1[2]<0) w1[2]=0;
    if (w1[3]<0) w1[3]=0;

    if (w2[0]<0) w2[0]=0;
    if (w2[1]<0) w2[1]=0;
    if (w2[2]<0) w2[2]=0;
    if (w2[3]<0) w2[3]=0;
    */

    //bug searching:  if the total of the weights is much greater than 1, exit...
    //assert(w1[0]+w1[1]+w1[2]+w1[3] < 1.1);
    //assert(w2[0]+w2[1]+w2[2]+w2[3] < 1.1);

    //diagnostics:
    if (tw1*tdiff > maxdt) {
      maxdt=tw1*tdiff;
    } else if (tw1*tdiff < mindt) {
      mindt=tw1*tdiff;	   
    }
    if (tw2*tdiff > maxdt) {
      maxdt=tw2*tdiff;
    } else if (tw2*tdiff < mindt) {
      mindt=tw2*tdiff;
    }
    if (tdiff > maxtdiff) {
      maxtdiff=tdiff;
    } else if (tdiff < mintdiff) {
      mintdiff=tdiff;
    }
    dtave+=tdiff;

#ifdef ASCII_OUT
    printf("(%8.2f, %8.2f), (", coords[i][0], coords[i][1]);
    for (long j=0; j<nchan-1; j++) printf("%8.2f, ", bt_interp[i][j]);
    printf("%8.2f)\n", bt_interp[i][nchan-1]);
#endif
  }

  //print the diagnostics:
  dtave/=n;
  fprintf(stderr, "diagnostic parameter          %8s   %8s   %8s\n",
	                     "min", "max", "average");
  fprintf(stderr, "time to scan line:           %10.3g %10.3g %10.3g\n",
		  mindt, maxdt, dtave/2);
  fprintf(stderr, "time interval:               %10.3g %10.3g %10.3g\n\n",
		  mintdiff, maxtdiff, dtave);

  //write the results to a file:
  fs=fopen(outfile, "w");
  fwrite(&nchan, sizeof(long), 1, fs);
  for (long i=0; i<n; i++) {
    fwrite(bt_interp[i], sizeof(float), nchan, fs);
  }
    
  fclose(fs);
  for (long i=0; i<n; i++) delete [] bt_interp[i];
  delete bt_interp;

  delete [] coords[0];
  delete [] coords;

  if (basepath != NULL) delete [] basepath;

}

