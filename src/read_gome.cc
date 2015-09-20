#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <regex.h>

#include "read_gome.h"
#include "time_class.h"
#include "peteys_tmpl_lib.h"

FILE *rg_logfile=NULL;

//generates a single spectra in a geometric series:
float * default_spectra(float lambda0, float lambdan, int n) {
  float *spec;
  spec=new float[n+1];
  for (long i=0; i<=n; i++) {
    spec[i]=lambda0*pow(lambdan/lambda0, (float) i/(float) n);
  }
  return spec;
}
  
//generates a single spectra in a geometric series:
float * default_spectra_mid(float lambda0, float lambdan, int n) {
  float *spec;
  spec=new float[n];
  for (long i=0; i<n; i++) {
    spec[i]=lambda0*pow(lambdan/lambda0, (i+0.5)/(float) n);
  }
  return spec;
}
  
#define NHEAD 13

#define MAXLL 200

#define NPPERSCAN 3
//#define NCHAN NPMD+5

#define NBANDMAX 4
#define NRECBAND 5

#define NSPECRAWMAX NBAND*1024

time_class convert_date(char *date) {
  time_class t;
  int year;
  int mon;
  int day;
  int hour;
  int min;
  float sec;

  sscanf(date, "%d2", &day);
  sscanf(date+7, "%d4", &year);
  sscanf(date+12, "%d2", &hour);
  sscanf(date+15, "%d2", &min);
  sscanf(date+18, "%f", &sec);

  date[6]=0;

  if (strcmp(date+3, "JAN")==0) {
    mon=1;
  } else if (strcmp(date+3, "FEB")==0) {
    mon=2;
  } else if (strcmp(date+3, "MAR")==0) {
    mon=3;
  } else if (strcmp(date+3, "APR")==0) {
    mon=4;
  } else if (strcmp(date+3, "MAY")==0) {
    mon=5;
  } else if (strcmp(date+3, "JUN")==0) {
    mon=6;
  } else if (strcmp(date+3, "JUL")==0) {
    mon=7;
  } else if (strcmp(date+3, "AUG")==0) {
    mon=8;
  } else if (strcmp(date+3, "SEP")==0) {
    mon=9;
  } else if (strcmp(date+3, "OCT")==0) {
    mon=10;
  } else if (strcmp(date+3, "NOV")==0) {
    mon=11;
  } else if (strcmp(date+3, "DEC")==0) {
    mon=12;
  } else {
    fprintf(stderr, "Unrecognized month: %s\n", date+3);
  }

  t.init(year, mon, day, hour, min, sec);

  return t;

}

#define SUN 0
#define GROUND 1
#define DUMMY 1

int NBAND=4;
float glob_specmin[NBANDMAX]={285, 313, 405, 602};
float glob_specmax[NBANDMAX]={313, 405, 602, 790};
int glob_bind[NBANDMAX]={0, 1, 2, 3};

//1. Why don't I pass this as a parameter??
//2. Do the bands selected have to match the bands originally extracted?
void rg_select_band(int b1b, int b2b, int b3, int b4) {
  int bflag[NBANDMAX];
  int k;
//  float specmin[NBANDMAX]={285, 313, 405, 602};
//  float specmax[NBANDMAX]={313, 405, 602, 790};

  bflag[0]=b1b; bflag[1]=b2b; bflag[2]=b3; bflag[3]=b4;

  k=0;
  for (int i=0; i<NBANDMAX; i++) {
    if (bflag[i]) {
      //glob_specmin[k]=specmin[i];
      //glob_specmax[k]=specmax[i];
      glob_bind[k]=i;
      k++;
    }
  }
  NBAND=k;

  printf("%d bands selected\n", NBAND);
}

#define NMATCH 10
#define REG_CFLAGS REG_EXTENDED
#define REG_EFLAGS 0
regmatch_t *rg_regex_match;

//compiled regular expressions:
regex_t *rg_date;
regex_t *rg_solar_spectrum;
regex_t *rg_channel;
regex_t *rg_earthshine;
regex_t *rg_ground_pixel;
regex_t *rg_band;
regex_t *rg_pmd;

//compile the regular expressions and allocate the match arrays:
int rg_compile_regex() {
  int r1, r2, r3, r4, r5, r6, r7;
  rg_date=new regex_t;
  r1=regcomp(rg_date, "([0-9]{2}\-[A-Z]{3}\-[0-9]{4} +[0-9]{2}\:[0-9]{2}\: ?[0-9]{1,2}\.[0-9]+)", REG_CFLAGS);
  rg_solar_spectrum=new regex_t;
  r2=regcomp(rg_solar_spectrum, "(Solar +Spectrum) +([0-9]{2}\-[A-Z]{3}\-[0-9]{4} +[0-9]{2}\:[0-9]{2}\: ?[0-9]{1,2}\.[0-9]+)", REG_CFLAGS);
  rg_channel=new regex_t;
  r3=regcomp(rg_channel, "(Channel +[0-9]) +([0-9]{3}\.[0-9]+) +([0-9]{3}\.[0-9]+) +([0-9]+) +", REG_CFLAGS);
  rg_earthshine=new regex_t;
  r4=regcomp(rg_earthshine, "(Earthshine +Spectrum) +([0-9]{2}\:[0-9]{2}\: ?[0-9]{1,2}\.[0-9]+) +([0-9]{2}\:[0-9]{2}\: ?[0-9]{1,2}\.[0-9]+) +([0-9]+)", REG_CFLAGS);
  rg_ground_pixel=new regex_t;
  r5=regcomp(rg_ground_pixel, "(Ground +Pixel) +([0-9]+) +([0-9]+) +([0-9]+)", REG_CFLAGS);
  rg_band=new regex_t;
  r6=regcomp(rg_band, "(Band +[0-9][ABab]?) +([0-9]+\.[0-9]+) +([0-9]+\.[0-9]+) +([0-9]+\.[0-9]+) +([0-9]+)", REG_CFLAGS);
  rg_pmd=new regex_t;
  r7=regcomp(rg_pmd, "(PMD) +([0-9]+) +([0-9]+) +([0-9]+)", REG_CFLAGS);

  if (r1 != 0 || r2 != 0 || r3 != 0 || r4 != 0 || r5 != 0 || r6 != 0 || r7 != 0) {
    fprintf(stderr, "%d %d %d %d %d %d %d\n", r1, r2, r3, r4, r5, r6, r7);
    return 1;
  }
  rg_regex_match=new regmatch_t[NMATCH];

  return 0;
}

void rg_delete_regex() {
  regfree(rg_date);
  delete rg_date;
  regfree(rg_solar_spectrum);
  delete rg_solar_spectrum;
  regfree(rg_channel);
  delete rg_channel;
  regfree(rg_earthshine);
  delete rg_earthshine;
  regfree(rg_ground_pixel);
  delete rg_ground_pixel;
  regfree(rg_band);
  delete rg_band;

  delete [] rg_regex_match;
}

//gets a single band from a single pixel (or solar spectrum)
float ** rg_get_band(FILE *in, int &n, char tp, char dum=0) {
  char line[MAXLL];
  char band[2];
  float integration_time;
  float spec0;
  float specf;
  float **data;
  int rex;

  //printf(line);

  //should use regular expressions rather than fixed columns...
  if (tp == SUN) {
    //sscanf(line+9, "%f %f %d", &spec0, &specf, &n);
    do {
      if (fgets(line, MAXLL, in)==NULL) {
        fprintf(stderr, "read_gome: read to end of file searching for solar spectra\n");
        return NULL;
      }
      rex=regexec(rg_channel, line, NMATCH, rg_regex_match, REG_EFLAGS);
    } while (rex!=0);
    sscanf(line+rg_regex_match[2].rm_so, "%g", &spec0);
    sscanf(line+rg_regex_match[3].rm_so, "%g", &specf);
    sscanf(line+rg_regex_match[4].rm_so, "%d", &n);
    line[rg_regex_match[1].rm_eo+1]='\0';
    if (rg_logfile!=NULL) fprintf(rg_logfile, "Reading %s: %d frequencies between %g and %g\n", line+rg_regex_match[1].rm_so, n, spec0, specf);
    //printf(line+9);
  } else {
    //sscanf(line+7, "%f %f %f %d", &integration_time, &spec0, &specf, &n);
    //printf("%f %f %f\n", integration_time, spec0, specf);
    do {
      if (fgets(line, MAXLL, in)==NULL) {
        fprintf(stderr, "read_gome: read to end of file searching for ground spectra\n");
        return NULL;
      }
      rex=regexec(rg_band, line, NMATCH, rg_regex_match, REG_EFLAGS);
    } while (rex!=0);
    
    sscanf(line+rg_regex_match[2].rm_so, "%g", &integration_time);
    sscanf(line+rg_regex_match[3].rm_so, "%g", &spec0);
    sscanf(line+rg_regex_match[4].rm_so, "%g", &specf);
    sscanf(line+rg_regex_match[5].rm_so, "%d", &n);
    line[rg_regex_match[1].rm_eo+1]='\0';
    if (rg_logfile!=NULL) fprintf(rg_logfile, "Reading %s: %d frequencies between %g and %g\n", line+rg_regex_match[1].rm_so, n, spec0, specf);
    //printf(line+9);
  }
  //printf("%d frequencies found in %6s\n", n, line);

  if (dum != 0) {
    for (long i=0; i<n; i++) fgets(line, MAXLL, in);
    return NULL;
  }

  data=new float *[NRECBAND];
  data[0]=new float[n*NRECBAND];
  for (int i=1; i<NRECBAND; i++) data[i]=data[0]+i*n;

  for (long i=0; i<n; i++) {
    fgets(line, MAXLL, in);
    sscanf(line, "%f %f %f %f %f", data[0]+i, data[1]+i, data[2]+i, data[3]+i, data[4]+i);
    //printf("%f\n", data[0]+i);
    //printf("%s", line);
  }

  return data;

}

//gets all bands and lumps them into a single spectrum in two arrays
int rg_get_raw(FILE *in, float *specraw, float *countsraw, char tp) {
  //we use four bands and only take them within these limits:
  //(we try to form a continuous spectra}

  int sind, find;	//bounds on the current band

  int nraw_ttl;
  float **specdata;
  int nraw;
  int nband2;
  int *bind;


  if (tp == SUN) {
    nband2=NBANDMAX;
    bind=new int[NBANDMAX];
    for (int i=0; i<NBANDMAX; i++) bind[i]=i;
  } else {
    nband2=NBAND;
    bind=glob_bind;
  }

  nraw_ttl=0;
  for (long j=0; j<nband2; j++) {
    specdata=rg_get_band(in, nraw, tp);
    if (specdata==NULL) return 0;		//empty pixel
    sind=bin_search(specdata[0], nraw, glob_specmin[bind[j]], -1);
    find=bin_search(specdata[0], nraw, glob_specmax[bind[j]], -1);
    for (int h=sind; h<=find; h++) {
      specraw[nraw_ttl]=specdata[0][h];
      //printf("%f\n", specraw[nraw_ttl]);
      countsraw[nraw_ttl]=specdata[1][h];
      nraw_ttl++;
    }
    delete [] specdata[0];
    delete [] specdata;
  }

  if (tp == SUN) delete [] bind;

  return nraw_ttl;

}

//interpolates the raw spectra into the desired spectra for whatever
//retrieval project...
void rg_interpolate_raw(float *specraw, float *countsraw, int nraw, float *spec, float *counts, int nchan) {
  double l1, l2;		//for interpolating the spectra
  long lastind;

  lastind=-1;

  //interpolate the raw sunshine spectra to our grid:
  l1=interpolate(specraw, nraw, spec[0], lastind);
  for (long j=0; j<nchan; j++) {
    l2=interpolate(specraw, nraw, spec[j+1], lastind);
    //printf("%f %f %d %d\n", spec[j], spec[j+1], l1, l2);
    counts[j]=(1+long(l1)-l1)*counts[long(l1)]+(l2-long(l2))*counts[long(l2)];
    for (long h=ceil(l1); h<floor(l2); h++) counts[j]+=countsraw[h];
    counts[j]/=(l2-l1);
    l1=l2;
  }
}

//reads GOME data from an ASCII file read from the level 1b data
//and interpolates it into our desired output spectra
//returns the data in a single structure
gome_data * read_gome(FILE *in, 		//ASCII input file
			float **spec,		//desired output spectra
			int * nchan,		//number of channels in each band
			int nb,			//number of bands
			int get_pmd) {		//PMD flag

  int npix;
  int ipi;
  int isc;
  int k;

  //FILE *in;
  char line[MAXLL];
  double vald[NPPERSCAN];
  float val[10];

  char datestr[32];

  gome_data *data;

  double missing=9.99999e-99;

  float lon;

  float *sunshine;	//sunshine spectra

  float countsraw[NSPECRAWMAX];
  float specraw[NSPECRAWMAX];
  int nraw;
  int ncum;
  int nspec;		//number of spectra, not including PMD

  time_class tsolar;
  char tstr[30];

  data=new gome_data;

  //in=fopen(fname, "r");

  if (rg_compile_regex() != 0) {
    fprintf(stderr, "read_gome: fatal error; failed to compile regular expressions\n");
    exit(1);
  }
 
  nspec=0;
  for (int i=0; i<nb; i++) nspec+=nchan[i];

  data->npix=0;
  if (get_pmd==1) data->nchan=nspec+NAUX+NPMD; else data->nchan=nspec+NAUX;
  data->missing=0./0.;

  //read in the sunshine spectra:
  //printf("Reading sun spectra\n");

  do {
    if (fgets(line, MAXLL, in)==NULL) {
      fprintf(stderr, "read_gome: failed to find solar spectra\n");
      rg_delete_regex();
      delete data;
      return NULL;
    }
    //printf("%s", line);
  } while (regexec(rg_solar_spectrum, line, NMATCH, rg_regex_match, REG_EFLAGS) != 0);

  tsolar=convert_date(line+rg_regex_match[2].rm_so);
  tsolar.write_string(tstr);
  line[rg_regex_match[1].rm_eo+1]='\0';
  if (rg_logfile!=NULL) fprintf(rg_logfile, "Reading %s at %s\n", line+rg_regex_match[1].rm_so, tstr);

  sunshine=new float[nspec];

  ncum=0;
  if (nb > 0) {
    nraw=rg_get_raw(in, specraw, countsraw, SUN);
    for (int j=0; j<nb; j++) {
      rg_interpolate_raw(specraw, countsraw, nraw, spec[j], sunshine+ncum, nchan[j]);
      ncum+=nchan[j];
    }
  } else {
    //do a dummy read:
    for (int j=0; j<NBANDMAX; j++) {
      rg_get_band(in, nraw, SUN, DUMMY);
    }
  }

  fgets(line, MAXLL, in);
  //sscanf(line+46, "%d", &npix);
  if (regexec(rg_earthshine, line, NMATCH, rg_regex_match, REG_EFLAGS) != 0) {
    fprintf(stderr, "read_gome: failed to match Earthshine header\n%s", line);
    rg_delete_regex();
    delete data;
    return NULL;
  }
  sscanf(line+rg_regex_match[4].rm_so, "%d", &npix);
  line[rg_regex_match[1].rm_eo+1]='\0';
  if (rg_logfile) fprintf(rg_logfile, "%s: %d ground pixels\n", line+rg_regex_match[1].rm_so, npix);

  data->data=new gome_rec[npix];

  for (long i=0; i<npix; i++) {

    if (fgets(line, MAXLL, in)==NULL) break;
    //sscanf(line+12, "%d %d %d", &ipi, &k, &isc);
    //printf(line);
    if (regexec(rg_ground_pixel, line, NMATCH, rg_regex_match, REG_EFLAGS) != 0) {
      fprintf(stderr, "read_gome: failed to match Ground Pixel header\n%s", line);
      continue;
    }
    sscanf(line+rg_regex_match[2].rm_so, "%d", &ipi);
    sscanf(line+rg_regex_match[3].rm_so, "%d", &k);
    sscanf(line+rg_regex_match[4].rm_so, "%d", &isc);
    line[rg_regex_match[1].rm_eo+1]='\0';
    if (rg_logfile!=NULL) fprintf(rg_logfile, "Reading %s: %d %d %d ", line+rg_regex_match[1].rm_so, i, ipi, isc);
    
    fgets(line, MAXLL, in);
    if (regexec(rg_date, line, NMATCH, rg_regex_match, REG_EFLAGS) != 0) {
      fprintf(stderr, "read_gome: failed to match date\n%s", line);
      break;
    }
    //printf("%s\n", line);
    data->data[data->npix].date=convert_date(line);
    data->data[data->npix].date.write_string(datestr);
    if (rg_logfile!=NULL) fprintf(rg_logfile, "measured on %s ", datestr);
    data->data[data->npix].counts=new float[data->nchan];

    fgets(line, MAXLL, in);
    fgets(line, MAXLL, in);

    fgets(line, MAXLL, in);
    sscanf(line, "%f %f %f %f", val, val+1, 
		    data->data[data->npix].counts,
		    data->data[data->npix].counts+1);

    fgets(line, MAXLL, in);
    sscanf(line, "%f %f %f %f", val, val+1, 
		    data->data[data->npix].counts+2,
		    data->data[data->npix].counts+3);

    fgets(line, MAXLL, in);
    fgets(line, MAXLL, in);
    fgets(line, MAXLL, in);
    fgets(line, MAXLL, in);
    //printf(line);
    sscanf(line, "%f %f %f %f %f %f %f %f %f %f", 
		    val, val+1, val+2, val+3, val+4, val+5, val+6, val+7, 
		    &data->data[data->npix].lat, &lon);
    if (lon > 180.) lon=lon-360.;
    data->data[data->npix].lon=lon;
    if (rg_logfile!=NULL) fprintf(rg_logfile, "at (%8.3f, %7.3f)\n", data->data[data->npix].lon, 
			data->data[data->npix].lat);

    //store the pixel number:
    data->data[data->npix].counts[4]=ipi;
    data->data[data->npix].counts[5]=isc;

    //fgets(line, MAXLL, in);
    //printf(line);
    
    //printf("(%f, %f) %s; %ld\n", data.data[data.nscan].lon[isc], 
	//	    data.data[data.nscan].lat[isc], datestr, isc);

    //read in and throw away (or not) the PMD data:
    if (get_pmd==1) {
      do {
        fgets(line, MAXLL, in);
        if (regexec(rg_band, line, 0, NULL, REG_EFLAGS) == 0) {
          fprintf(stderr, "read_gome: failed to find PMD data\n");
          fseek(in, -strlen(line), SEEK_CUR);
          get_pmd=0;
          break;
        }
      } while (regexec(rg_pmd, line, NMATCH, rg_regex_match, REG_EFLAGS)!=0);
      if (get_pmd==0) for (long j=0; j<16; j++) {
        fgets(line, MAXLL, in);
        fscanf(in, "%lg%lg%lg", vald, vald+1, vald+2);
        if (vald[0] == missing) data->data[data->npix].counts[NAUX+j*3]=data->missing;
                else data->data[data->npix].counts[NAUX+j*3]=vald[0];
        if (vald[1] == missing) data->data[data->npix].counts[NAUX+j*3+1]=data->missing;
                else data->data[data->npix].counts[NAUX+j*3+1]=vald[1];
        if (vald[2] == missing) data->data[data->npix].counts[NAUX+j*3+2]=data->missing;
                else data->data[data->npix].counts[NAUX+j*3+2]=vald[2];
        //printf(line);
      } 
    } //else {
      //for (long j=0; j<16; j++) fgets(line, MAXLL, in);
    //}

    if (nb > 0) {
      int offset;
      if (get_pmd==1) offset=NAUX+NPMD; else offset=NAUX;
      nraw=rg_get_raw(in, specraw, countsraw, GROUND);
      if (nraw==0) continue;		//empty pixel
      ncum=offset;
      for (int j=0; j<nb; j++) {
        rg_interpolate_raw(specraw, countsraw, nraw, spec[j], data->data[data->npix].counts+ncum, nchan[j]);
        ncum+=nchan[j];
      }
      //divide by the sunshine spectra:
      for (int j=0; j<nspec; j++) data->data[data->npix].counts[j+offset]/=sunshine[j];
    } else {
      //do a dummy read:
      for (int j=0; j<NBAND; j++) {
        rg_get_band(in, nraw, GROUND, DUMMY);
        if (nraw==0) break;
      }
      if (nraw==0) continue;		//empty pixel
    }

    data->npix++;

  }

  printf("%d pixels extracted\n", data->npix);

  delete [] sunshine;
  rg_delete_regex();

  return data;


}

void delete_gome_data(gome_data *data) {
  for (int4_t i=0; i<data->npix; i++) {
    delete [] data->data[i].counts;
  }

  delete [] data->data;
  delete data;

}

#define COMMENT_LENGTH 500

//***note: the 'swap_end' parameter exists for compatibility reasons only
//	   and currently doesn't do anything

int4_t rawread_gome_head(char *filename, time_class & t1, time_class & t2,
		int4_t &nchan, int swap_end) {
  FILE *fs;
  char comment[COMMENT_LENGTH];
  int4_t nscan;

  int4_t year;
  double day;

  //open the file:
  fs=fopen(filename, "r");
  if (fs == NULL) {
    fprintf(stderr, "File, %s, could not be opened for reading\n", filename);
    nchan=0;
    return -1;
  }

  //read the comment header:
  fread(comment, sizeof(char), COMMENT_LENGTH, fs);
  fprintf(stderr, "%s", comment);

  //read the attributes:
  fread(&year, sizeof(year), 1, fs);
  fread(&day, sizeof(day), 1, fs);
  //this is awkward:
  t1.init(year, 1, 1, 0, 0, 0);
  t1.add(day);
 
  fread(&year, sizeof(year), 1, fs);
  fread(&day, sizeof(day), 1, fs);
  t2.init(year, 1, 1, 0, 0, 0);
  t2.add(day);

  fread(&nscan, sizeof(nscan), 1, fs);
  fread(&nchan, sizeof(nchan), 1, fs);
  //fread(&missing, sizeof(missing), 1, fs);
 
  fclose(fs);

  return nscan;

}

int4_t rawwrite_gome(char * filename, gome_data *data, int swap_end) {
  FILE *fs;
  int4_t nwrit;
  char comment[COMMENT_LENGTH];
  char tstring1[30], tstring2[30];

  int4_t year;
  double day;

  //open the file:
  fs=fopen(filename, "w");
  if (fs == NULL) {
    fprintf(stderr, "File, %s, could not be opened for writing\n", filename);
    return 1;
  }

  //write the headers:
 
  //comment header:
  data->data[0].date.write_string(tstring1);
  data->data[data->npix-1].date.write_string(tstring2);
  sprintf(comment, "GOME or GOME-derived data in raw binary format\nDate range:%30s-%30s\n%10d scan lines; %3d channels\n", 
		tstring1, tstring2, data->npix, data->nchan);
  fwrite(comment, sizeof(char), COMMENT_LENGTH, fs);

  //attributes:
  year=(int4_t) data->data[0].date.year();
  day=data->data[0].date.doy();
  fwrite(&year, sizeof(year), 1, fs);
  fwrite(&day, sizeof(day), 1, fs);

  year=(int4_t) data->data[data->npix-1].date.year();
  day=data->data[data->npix-1].date.doy();
  fwrite(&year, sizeof(year), 1, fs);
  fwrite(&day, sizeof(day), 1, fs);

  fwrite(&data->npix, sizeof(data->npix), 1, fs);
  fwrite(&data->nchan, sizeof(data->nchan), 1, fs);
  fwrite(&data->missing, sizeof(data->missing), 1, fs);

  //write out the data scan-line by scan-line:
  for (int4_t i=0; i<data->npix; i++) {
    year=(int4_t) data->data[i].date.year();
    day=data->data[i].date.doy();
    //try to align everything with the word boundaries
    //for easy byte-swapping:
    fwrite(&year, sizeof(year), 1, fs);
    fwrite(&day, sizeof(day), 1, fs);

    fwrite(&data->data[i].lon, sizeof(float), 1, fs);
    fwrite(&data->data[i].lat, sizeof(float), 1, fs);
    fwrite(data->data[i].counts, sizeof(float), data->nchan, fs);
  }

  fclose(fs);

  return 0;

}
     
gome_data * rawread_gome(char * filename, int swap_end) {
  FILE *fs;
  time_class t1, t2;
  char comment[COMMENT_LENGTH];
  gome_data *data;

  int4_t year;
  double day;

  //open the file:
  fs=fopen(filename, "r");
  if (fs == NULL) {
    fprintf(stderr, "File, %s, could not be opened for reading\n", filename);
    return NULL;
  }
	  
  //read the comment header:
  fread(comment, sizeof(char), COMMENT_LENGTH, fs);
  fprintf(stderr, "%s", comment);

  data=new gome_data;

  //read the attributes:
  fread(&year, sizeof(year), 1, fs);
  fread(&day, sizeof(day), 1, fs);

  fread(&year, sizeof(year), 1, fs);
  fread(&day, sizeof(day), 1, fs);

  fread(&data->npix, sizeof(data->npix), 1, fs);
  fread(&data->nchan, sizeof(data->nchan), 1, fs);
  fread(&data->missing, sizeof(data->missing), 1, fs);

  data->data=new gome_rec[data->npix];

  //read the data pixel by pixel:
  for (long i=0; i<data->npix; i++) {
    //allocate space:
    fread(&year, sizeof(year), 1, fs);
    fread(&day, sizeof(day), 1, fs);
    data->data[i].date.init(year, 1, 1, 0, 0, 0);
    data->data[i].date.add(day);

    fread(&data->data[i].lon, sizeof(float), 1, fs);
    fread(&data->data[i].lat, sizeof(float), 1, fs);

    data->data[i].counts=new float[data->nchan];
    fread(data->data[i].counts, sizeof(float), data->nchan, fs);
  }

  fclose(fs);

  return data;

}
