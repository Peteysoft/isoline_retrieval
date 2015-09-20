#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "peteys_tmpl_lib.h"
#include "colloc_obj.h"

#include "global_metric.h"

#define MAX_LL 200
#define SWAP_ENDIAN 1

#define COMMENT_LENGTH 200

//#define CHECK_INDEX_FILE 

colloc_obj::colloc_obj() {
  amsu_head_reader=&rawread_gome_head;
  amsu_reader=&rawread_gome;

  filelist=NULL;
  amsu_data=NULL;
  //tstart=NULL;
  //tend=NULL;
  nfiles=0;
  fname=NULL;
}

colloc_obj::colloc_obj(char *initfile, long mnl, char *basepath1) {
  char *filename;
  FILE *fs;

  char line[MAX_LL];
  long ll;			//line length
  char ** flist1;		//for sorting file-list
  time_class *ts;
  time_class *t1;
  long *sind;			//sorting indices

  //throw-away:
  int4_t nchan;
  int4_t *npix;
  long nlen;		//length of filename
  long baselen;		//length of base path name

  char tstring[30], ts1[30], ts2[30];

  amsu_head_reader=&rawread_gome_head;
  amsu_reader=&rawread_gome;

  if (basepath1 == NULL) {
    baselen=0;
    basepath=NULL;
    strcpy(filename, "");
  } else {
    baselen=strlen(basepath1);
    basepath=new char[baselen+1];
    strcpy(basepath, basepath1);
    strcpy(filename, basepath1);
  }

  fs=fopen(initfile, "r");

  //mark our position:
  //count the lines:
  for (nfiles=0; feof(fs) == 0; nfiles++) {
    if (fgets(line, MAX_LL, fs)==NULL) break;
    if (strlen(line) <= 1) break;
  }

  printf("%d files found in initialisation file %s\n", nfiles, initfile);

  //rewind the file and read the files in earnest:
  fseek(fs, 0, SEEK_SET);
  filelist=new char *[nfiles];
  tstart=new time_class[nfiles];
  tend=new time_class[nfiles];
  npix=new int4_t[nfiles];
  for (long i=0; i<nfiles; i++) {
    //read the name of the file:
    fgets(line, MAX_LL, fs);
    ll=strlen(line);
    line[ll-1]=0;
    filelist[i]=new char[ll];
    strcpy(filelist[i], line);

    //read start and end dates:
    strcpy(filename+baselen, filelist[i]);
    npix[i]=(*amsu_head_reader)(filename, tstart[i], tend[i], nchan, swap_end_flag);
    ntotal+=npix[i];
  }

  fclose(fs);

  //sort the file list and dates by start date:
  sind=heapsort(tstart, nfiles);
  t1=map_vector(tstart, sind, nfiles);
  delete [] tstart;
  tstart=t1;

  t1=map_vector(tend, sind, nfiles);
  delete [] tend;
  tend=t1;

  flist1=new char *[nfiles];
  for (long i=0; i<nfiles; i++) flist1[i]=filelist[sind[i]];
  delete [] filelist;
  filelist=flist1;

  delete [] sind;

  printf("The following files' start and end dates were loaded:\n");
  for (long i=0; i<nfiles; i++) {
    tstart[i].write_string(ts1);
    tend[i].write_string(ts2);
    printf("%s: %s-%s\n", filelist[i], ts1, ts2);
  }

  //set number of loaded files and max. allowed:
  nloaded=0;
  maxload=mnl;

  amsu_data=new gome_data *[nfiles];
  for (long i=0; i<nfiles; i++) amsu_data[i]=NULL;

  tind=NULL;
  t=NULL;
  nind=-1;
}

colloc_obj::~colloc_obj() {
//  printf("In: ~amsu_interp_obj\n");
  for (long i=0; i<nfiles; i++) delete [] filelist[i];
  delete [] filelist;

  //delete [] tstart;
  //delete [] tend;

  for (long i=0; i<nfiles; i++) if (amsu_data[i] != NULL) delete_gome_data(amsu_data[i]);
  delete [] amsu_data;

}

//index the data by date:
int colloc_obj::index_data() {
  char *filename;
  int *npix;
  int k;
  gome_data data;
  double tref_d;

  float *t1;
  short *tind1;
  long *sind;
  int ndup;

  if (basepath == NULL) {
    baselen=0;
    strcpy(filename, "");
  } else {
    baselen=strlen(basepath);
    strcpy(filename, basepath);
  }

  nind=0;
  for (int i=0; i<nfiles; i++) {
    strcpy(filename+baselen, filelist[i]);
    npix=(*amsu_head_reader)(filename, tstart[i], tend[i], nchan, swap_end_flag);
    nind+=npix;
  }

  tind1=new short[nind*2];
  t1=new float[nind];
  tref=tstart[0];
  tref_d=(double) tref;

  k=0;
  for (short i=0; i<nfiles; i++) {
    tind1[k*2]=i;
    strcpy(filename+baselen, filelist[i]);
    (*amsu_reader) (filename, &data, swap_end_flag);
    for (short j=0; j<data->npix; j++) {
      t1[k]=(double) data->data[j].time-tref_d;
      tind1[k*2+1]=j;
    }
    k++;
  }

  //sort the dates, removing duplicates:
  sind=heapsort(t1, nind);
  if (t!=NULL) delete [] t;
  if (tind!=NULL) delete [] tind;
  t=new float[nind];
  tind=new short[nind*2];
  ndup=0;
  t[0]=t1[sind[0]];
  tind[0]=tind1[sind[0]*2];
  tind[1]=tind1[sind[0]*2+1];
  for (int i=1; i<nind; i++) {
    t[i-ndup]=t1[sind[i]];
    if (t[i]=t[i-1]) {
      ndup++;
      continue;
    }
    tind[(i-ndup)*2]=tind1[sind[i]*2];
    tind[(i-ndup)*2+1]=tind1[sind[i]*2+1];
  }
  delete [] t1;
  delete [] tind1;
  delete [] sind;
  delete_gome_data(&data);

  nind=nind-ndup;

  return 0;

}

int colloc_obj::write_indexfile(char *filename) {
  FILE *fs;
  char header[COMMENT_LENGTH];
  int32_t slen;
  int32_t amsu_read_type=RAW;

  fs=fopen(filename);

  sprintf(header, "Indexes dates for nadir-looking satellite data\n");
  fwrite(fs, sizeof(char), COMMENT_LENGTH, header);
  slen=strlen(basepath)+1;
  fwrite(fs, sizeof(slen), 1, &slen);
  fwrite(fs, sizeof(char), slen, basepath);
  fwrite(fs, sizeof(amsu_read_type), 1, &amsu_read_type);
  fwrite(fs, sizeof(nfiles), 1, &nfiles);
  for (int i=0; i<nfiles; i++) {
    slen=strlen(filelist[i])+1;
    fwrite(fs, sizeof(slen), 1, &slen);
    fwrite(fs, sizeof(char), slen, filelist[i]);
  }
  fwrite(fs, sizeof(tref), 1, &tref);
  fwrite(fs, sizeof(nind), 1, &nind);
  fwrite(fs, sizeof(float), nind, t);
  fwrite(fs, sizeof(short), nind*2, tind);
}

int colloc_obj::read_flist(char *initfile) {
  char *filename;
  FILE *fs;

  char line[MAX_LL];
  long ll;			//line length
  char ** flist1;		//for sorting file-list
  time_class *ts;
  time_class *t1;
  long *sind;			//sorting indices

  //throw-away:
  int4_t nchan;
  int4_t *npix;
  long nlen;		//length of filename
  long baselen;		//length of base path name

  char tstring[30], ts1[30], ts2[30];

  if (basepath == NULL) {
    baselen=0;
    strcpy(filename, "");
  } else {
    baselen=strlen(basepath);
    strcpy(filename, basepath1);
  }

  fs=fopen(initfile, "r");

  //mark our position:
  //count the lines:
  for (nfiles=0; feof(fs) == 0; nfiles++) {
    if (fgets(line, MAX_LL, fs)==NULL) break;
    if (strlen(line) <= 1) break;
  }

  printf("%d files found in initialisation file %s\n", nfiles, initfile);

  //rewind the file and read the files in earnest:
  fseek(fs, 0, SEEK_SET);
  filelist=new char *[nfiles];
  tstart=new time_class[nfiles];
  tend=new time_class[nfiles];
  npix=new int4_t[nfiles];
  for (long i=0; i<nfiles; i++) {
    //read the name of the file:
    fgets(line, MAX_LL, fs);
    ll=strlen(line);
    line[ll-1]=0;
    filelist[i]=new char[ll];
    strcpy(filelist[i], line);

    //read start and end dates:
    strcpy(filename+baselen, filelist[i]);
    npix[i]=(*amsu_head_reader)(filename, tstart[i], tend[i], nchan, swap_end_flag);
    ntotal+=npix[i];
  }

  fclose(fs);

  //sort the file list and dates by start date:
  sind=heapsort(tstart, nfiles);
  t1=map_vector(tstart, sind, nfiles);
  delete [] tstart;
  tstart=t1;

  t1=map_vector(tend, sind, nfiles);
  delete [] tend;
  tend=t1;

  flist1=new char *[nfiles];
  for (long i=0; i<nfiles; i++) flist1[i]=filelist[sind[i]];
  delete [] filelist;
  filelist=flist1;

  delete [] sind;

  printf("The following files' start and end dates were loaded:\n");
  for (long i=0; i<nfiles; i++) {
    tstart[i].write_string(ts1);
    tend[i].write_string(ts2);
    printf("%s: %s-%s\n", filelist[i], ts1, ts2);
  }
}

float colloc_obj::colloc_forward(float lon0,           //lon-lat coords of int. point
                        float lat0,
                        time_class t0,                  //time of interpolation point
                        float critdist,                 //maximum distance to satellite pixel
                        time_class &tforward,            //time of satellite measurement
			long &pixnum,
			char *&filename,
			long &nchan,
                        float *&val)                   //four data values
{

  //interpolation point:
  int ind0;
  short indf, inds;

  //amsu data:
  gome_data *data;

  //distances from each point in the scan line:
  float d;

  float t0f;	//test time as floating point

  //for output:
  char tstring[30];

  int ffound=-1;   //found forwards, found backwards

  if (t==NULL) index_data();

  t0f=(double) t0 - (double) tref;
  ind0=bin_search(t, nind, t0f, -1);

  //fs=fopen("testint.txt", "w");
  for (long i=ind0; i<nfiles; i++) {
    indf=tind[ind*2];
    inds=tind[ind*2+1];
    data=get_data_forward(indf);

    d=sdist(data->data[inds].lon, data->data[inds].lat, lon0, lat0);

    //we may have found our point 
    if (d < critdist) {
      //printf("dmin=%f; critdist=%f\n", dmin, critdist);
      pixnum=inds;
      filename=filelist[indf];
      tforward=data->data[inds].date;
      //this is very bad (but shouldn't affect our results...):
      val=data->data[inds].counts;
      ffound=d;
      break;
    }
  }

#ifdef PRINT_FLOC
  printf("f: %s %d %d\n", filelist[indf], inds, minind);
#endif

  nchan=data->nchan;

  //fclose(fs);

  return ffound;

}

float colloc_obj::colloc_backward(float lon0,           //lon-lat coords of int. point
                        float lat0,
                        time_class t0,                  //time of interpolation point
                        float critdist,                 //maximum distance to satellite pixel
                        time_class &tbackward,            //time of satellite measurement
			long &pixnum,
			char *& filename,
			long &nchan,
                        float *& val)                   //four data values
{

  //interpolation point:
  int ind0;
  short indf, inds;

  //amsu data:
  gome_data *data;

  //distances from each point in the scan line:
  float d;

  //for output:
  char tstring[30];

  int ffound=-1;   //found forwards, found backwards

  if (t==NULL) index_data();

  t0f=(double) t0 - (double) tref;
  ind0=bin_search(t, nind, t0f, -1);

  for (long i=ind0; i>=nind; i--) {
    indf=tind[i*2];
    inds=tind[i*2+1];
    data=get_data_backward(indf);

    d=sdist(data->data[inds].lon, data->data[inds].lat, lon0, lat0);

    //we may have found our point, 
    if (d < critdist) {
      //printf("dmin=%f; critdist=%f\n", dmin, critdist);
      pixnum=inds;
      filename=filelist[indf];
      tbackward=data->data[inds].date;
      val=data->data[inds].counts;
      ffound=d;
      break;
    }
  }

#ifdef PRINT_FLOC
  printf("f: %s %d %d\n", filelist[indf], inds, minind);
#endif

  nchan=data->nchan;

  //fclose(fs);

  return ffound;

}

