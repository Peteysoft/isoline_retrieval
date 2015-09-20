#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "peteys_tmpl_lib.h"
#include "amsu_interp_obj.h"

#include "global_metric.h"

#define MAX_LL 200
#define SWAP_ENDIAN 1

#define COMMENT_LENGTH 200

//#define CHECK_INDEX_FILE 

amsu_interp_obj::amsu_interp_obj() {
  filelist=NULL;
  amsu_data=NULL;
  //tstart=NULL;
  //tend=NULL;
  nfiles=0;
}

amsu_interp_obj::amsu_interp_obj(char *initfile, long mnl, char *basepath) {

  FILE *fs;
  FILE *fscheck;
  char checkfile[MAX_LL];

  char line[MAX_LL];
  long fp_amsustart;		//start of list of amsu files
  long ll;			//line length
  char ** flist1;		//for sorting file-list
  time_class *ts;
  long *sind;			//sorting indices

  //throw-away:
  long np, nchan;
  char comment[COMMENT_LENGTH];
  long nlen;		//length of filename
  long baselen;		//length of base path name
  long lonind, latind;
  long nrec;
  float lon, lat;
  float check_lon, check_lat;

  char tstring[30];

  int amsu_io_type;

  fs=fopen(initfile, "r");

  fread(comment, sizeof(char), COMMENT_LENGTH-sizeof(float), fs);
  printf("%s\n", comment);
  fread(&lon_offset, sizeof(float), 1, fs);
  lon_offset=lon_offset+180;
  fread(&amsu_io_type, sizeof(amsu_io_type), 1, fs);
  fread(&swap_end_flag, sizeof(swap_end_flag), 1, fs);

  //read the base file path:
  nlen=1;
  while ((char) fgetc(fs) != '\0') nlen++;
  if (basepath==NULL) {
    basepath=new char[nlen];
    baselen=nlen-1;
    fseek(fs, -nlen, SEEK_CUR);
    fread(basepath, sizeof(char), nlen, fs);
  } else {
    baselen=strlen(basepath);
  }

  fread(&nfiles, sizeof(nfiles), 1, fs);
//#ifdef CHECK_INDEX_FILE
  if (amsu_io_type == RAW) {
    printf("%4d RAW files\n", nfiles);
  } else {
    printf("%4d NORMAL files\n", nfiles);
  }
  printf("%s\n", basepath);
//#endif

  filelist=new char *[nfiles];
  //printf("%d file names found in index file %s\n", nfiles, initfile);
  for (long i=0; i<nfiles; i++) {
    nlen=1;
    while ((char) fgetc(fs) != '\0') nlen++;
    filelist[i]=new char[nlen+baselen+1];
    strcpy(filelist[i], basepath);
    fseek(fs, -nlen, SEEK_CUR);
    fread(filelist[i]+baselen, sizeof(char), nlen, fs);
//#ifdef CHECK_INDEX_FILE
    printf("%s\n", filelist[i]+baselen);
//#endif
  }

  fread(&bin_nlat, sizeof(long), 1, fs);
  bin_nlon=new long[bin_nlat+2];
  fread(bin_nlon, sizeof(long), bin_nlat+2, fs);

  //printf("Reading bins: \n");
  bin=new amsu_index_bin *[bin_nlat+2];

  for (long i=0; i<bin_nlat+2; i++) {
    printf("%d\n", bin_nlon[i]);
    bin[i]=new amsu_index_bin[bin_nlon[i]];
	    
    for (long j=0; j<bin_nlon[i]; j++) {
	    
      check_lon=j*360./bin_nlon[i]-lon_offset;
      check_lat=i*180./bin_nlat-90;

      fread(&latind, sizeof(long), 1, fs);
      fread(&lonind, sizeof(long), 1, fs);
      fread(&lon, sizeof(float), 1, fs);
      fread(&lat, sizeof(float), 1, fs);
      fread(&nrec, sizeof(long), 1, fs);
      bin[i][j].n=nrec;

      if (latind != i || lonind != j) {
        fprintf(stderr, "File is corrupted: %d vs. %d, %d vs %d\n", latind, i, lonind, j);
	fprintf(stderr, "(%f, %f) vs. (%f, %f)\n", lon, lat, check_lon, check_lat);
	exit(26);
      }

#ifdef CHECK_INDEX_FILE
      sprintf(checkfile, "check_bins/check_bin(%7.2f,%7.2f).txt", lon, lat);
      fscheck=fopen(checkfile, "w");
#endif
      printf("(%6d, %6d); (%8.2f, %8.2f): %d\n", latind, lonind, lon, lat, nrec);
      //printf("%6d (%8.2f, %8.2f)\n", check_binind, check_lon, check_lat);

      if (nrec == 0) {
        continue;
      }
      bin[i][j].date=new time_class[nrec];
      bin[i][j].file_index=new long[nrec];
      bin[i][j].scan_index=new long[nrec];
      fread(bin[i][j].date, sizeof(time_class), nrec, fs);
      fread(bin[i][j].file_index, sizeof(long), nrec, fs);
      fread(bin[i][j].scan_index, sizeof(long), nrec, fs);

#ifdef CHECK_INDEX_FILE
      for (long k=0; k<nrec; k++) {
        bin[i][j].date[k].write_string(tstring);
        //printf("%29s %4d %10d\n", tstring, bin[i][j].file_index[k], bin[i][j].scan_index[k]);
        fprintf(fscheck, "%29s %4d %10d\n", tstring, bin[i][j].file_index[k], bin[i][j].scan_index[k]);
      }
#endif

#ifdef CHECK_INDEX_FILE
      fclose(fscheck);
#endif
    }

    //set limit for polar bins:
    if (i == 1) {
      polar_latthresh=90+lat;
      printf("threshold for polar bins: %f", polar_latthresh);
    }
	    
  }

  fclose(fs);

  if (amsu_io_type==RAW) {
    amsu_head_reader=&rawread_amsu_1c_head;
    amsu_reader=&rawread_amsu_1c;
  } else {
    amsu_head_reader=&read_amsu_1c_head;
    amsu_reader=&read_amsu_1c;
  }

  //set number of loaded files and max. allowed:
  nloaded=0;
  maxload=mnl;

  amsu_data=new amsu_1c_data *[nfiles];
  for (long i=0; i<nfiles; i++) amsu_data[i]=NULL;

  delete [] basepath;

}

amsu_interp_obj::~amsu_interp_obj() {
//  printf("In: ~amsu_interp_obj\n");
  for (long i=0; i<nfiles; i++) delete [] filelist[i];
  delete [] filelist;

  //delete [] tstart;
  //delete [] tend;

  for (long i=0; i<nfiles; i++) if (amsu_data[i] != NULL) delete_amsu_1c_data(amsu_data[i]);
  delete [] amsu_data;

  for (long i=0; i<bin_nlat; i++) {
    for (long j=0; j<bin_nlon[i]; j++) {
      delete [] bin[i][j].date;
      delete [] bin[i][j].file_index;
      delete [] bin[i][j].scan_index;
    }
    delete [] bin[i];
  }

  delete [] bin;

}

int amsu_interp_obj::binind(float lon, float lat, long &lonind, long &latind) {
  //calculate the bin index:
  latind=(long)((lat+90.-polar_latthresh)/(180.-2*polar_latthresh)*bin_nlat)+1;
  lonind=(long)((lon+lon_offset)/360*bin_nlon[latind]);
}

#define MINW_EXTRAP -1	

int amsu_interp_obj::int_forward(float lon0,           //lon-lat coords of int. point
                        float lat0,
                        time_class t0,                  //time of interpolation point
                        float critdist,                 //maximum distance to satellite pixel
                        time_class &tforward,            //time of satellite measurement
			long &nchan,
                        float *val[4],                   //four data values
                        float weight[4],              //their weights
			int no_extrap)
{

//  long np, nchan;

  FILE *fs, *fs2;

  //bin index:
  long lonind, latind;

  long nrec;		//number of records in the bin
  amsu_index_bin *cur_bin;

  //interpolation point:
  long ind, indf, inds;

  //for looking backwards:
  time_class tnext;

  //amsu data:
  amsu_1c_data *data;

  float minw;
  float minw_best;	//highest minimum weight
  float *testv[4];	//equivalent measurements
  float testw[4];

  //distances from each point in the scan line:
  float d;
  float dmin;		//minimum distance
  long minind;		//where in the scan line...

  //for output:
  char tstring[30];
  char tstring1[30];

  int ffound=1;   //found forwards, found backwards

  //calculate the bin index:
  binind(lon0, lat0, lonind, latind);

  cur_bin=&bin[latind][lonind];
  nrec=bin[latind][lonind].n;

  if (nrec == 0) {
    fprintf(stderr, "Error: index bin (%f, %f) is empty\n", lon0, lat0);
    return ffound;
  }
  if (t0 > cur_bin->date[nrec-1]) {
    t0.write_string(tstring);
    cur_bin->date[nrec-1].write_string(tstring1);  
    fprintf(stderr, "Error: date (%s) later than any in bin (%s)\n", tstring, tstring1);
    return ffound;
  }  

  ind=bin_search(cur_bin->date, nrec, t0, -1)+1;

  minw_best=MINW_EXTRAP;

  //fs=fopen("testint.txt", "w");
  for (long i=ind; i<nrec; i++) {
    indf=cur_bin->file_index[i];
    inds=cur_bin->scan_index[i];

    data=get_data_forward(indf);

    dmin=sdist(data->data[inds].lon[0], data->data[inds].lat[0], lon0, lat0);
    minind=0;
    for (long j=1; j<data->np; j++) {
      d=sdist(data->data[inds].lon[j], data->data[inds].lat[j], lon0, lat0);
      if (d<dmin) {
        dmin=d;
        minind=j;
      }
    }
    //we may have found our point: get the weights, 
    //if none of them are negative, set the flag and break out of the loop
    //since we don't want to do an extrapolation...
    if (dmin < critdist) {
      //printf("dmin=%f; critdist=%f\n", dmin, critdist);
      derive_amsu_int_weights(data, inds, minind, lon0, lat0, testv, testw);
      //don't want any missing data:
      //if (testv[0] > data->missing && testv[1] > data->missing && 
	//		testv[2] > data->missing && testv[3] > data->missing) {
        if ((testw[0] >= 0 && testw[1]>=0 && testw[2]>=0 && testw[3]>=0) || no_extrap == 0) {
          tforward=cur_bin->date[i];
          weight[0]=testw[0]; weight[1]=testw[1]; weight[2]=testw[2]; weight[3]=testw[3];
          val[0]=testv[0]; val[1]=testv[1]; val[2]=testv[2]; val[3]=testv[3];
          ffound=0;
          break;
        } else {
          //in the event we don't find a set of four points enclosing the test
          //point, search for the set that requires the least extrapolation:
          minw=testw[0];
          for (int k=1; k<4; k++) if (testw[k]<minw) minw=testw[k];
          if (minw > minw_best) {
            tforward=cur_bin->date[i];
            weight[0]=testw[0]; weight[1]=testw[1]; weight[2]=testw[2]; weight[3]=testw[3];
	    val[0]=testv[0]; val[1]=testv[1]; val[2]=testv[2]; val[3]=testv[3];
            minw_best=minw;
	    ffound=-1;
	  }
          //printf("Negative weights!\n");
	}
      //}
    }
    //date[i].write_string(tstring);
    //fprintf(fs, "%24s %10.3f %10.3f %3d\n", tstring, critdist, dmin, minind);
  }

#ifdef PRINT_FLOC
  printf("f: %s %d %d\n", filelist[indf], inds, minind);
#endif

  nchan=data->nchan;

  //fclose(fs);

  return ffound;

}

int amsu_interp_obj::int_backward(float lon0,           //lon-lat coords of int. point
                        float lat0,
                        time_class t0,                  //time of interpolation point
                        float critdist,                 //maximum distance to satellite pixel
                        time_class &tbackward,            //time of satellite measurement
			long &nchan,
                        float *val[4],                   //four data values
                        float weight[4],              //their weights
			int no_extrap)
{

//  long np, nchan;

  FILE *fs, *fs2;

  //bin index:
  long lonind, latind;
  //long binind;

  long nrec;		//number of records in the bin
  amsu_index_bin *cur_bin;

  //interpolation point:
  long ind, indf, inds;

  //for looking backwards:
  time_class tnext;

  //amsu data:
  amsu_1c_data *data;

  float minw;
  float minw_best;	//highest minimum weight
  float *testv[4];	//equivalent measurements
  float testw[4];

  //distances from each point in the scan line:
  float d;
  float dmin;		//minimum distance
  long minind;		//where in the scan line...

  //for output:
  char tstring[30];
  char tstring1[30];

  int ffound=1;   //found forwards, found backwards

  //calculate the bin index:
  binind(lon0, lat0, lonind, latind);

  cur_bin=&bin[latind][lonind];
  nrec=bin[latind][lonind].n;

  if (nrec == 0) {
    fprintf(stderr, "Error: index bin (%f, %f) is empty\n", lon0, lat0);
    return ffound;
  }
  if (t0 < cur_bin->date[0]) {
    t0.write_string(tstring);
    cur_bin->date[0].write_string(tstring1);  
    fprintf(stderr, "Error: date (%s) before any in bin (%s)\n", tstring, tstring1);
    return ffound;
  }  

  ind=bin_search(cur_bin->date, cur_bin->n, t0, -1);
  //printf("%d entries found in bin #%d.  Starting at %d\n", 
	//	  bin[binind].n, binind, ind);

  minw_best=MINW_EXTRAP;

  //fs=fopen("testint.txt", "w");
  for (long i=ind; i>=0; i--) {
    indf=cur_bin->file_index[i];
    inds=cur_bin->scan_index[i];

    data=get_data_backward(indf);

    dmin=sdist(data->data[inds].lon[0], data->data[inds].lat[0], lon0, lat0);
    minind=0;
    for (long j=1; j<data->np; j++) {
      d=sdist(data->data[inds].lon[j], data->data[inds].lat[j], lon0, lat0);
      if (d<dmin) {
        dmin=d;
        minind=j;
      }
    }

    //printf("%d: dmin=%f\n", i, dmin);
    //we found our point: get the weights, set the flag and break out of the loop:
    if (dmin < critdist) {
      //printf("dmin=%f; critdist=%f\n", dmin, critdist);
      derive_amsu_int_weights(data, inds, minind, lon0, lat0, testv, testw);
      //don't want any missing data:
      //if (testv[0] > data->missing && testv[1] > data->missing && 
	//		testv[2] > data->missing && testv[3] > data->missing) {
        if ((testw[0] >= 0 && testw[1]>= 0 && testw[2]>=0 && testw[3]>=0) || no_extrap == 0) {
          tbackward=cur_bin->date[i];
          weight[0]=testw[0]; weight[1]=testw[1]; weight[2]=testw[2]; weight[3]=testw[3];
          val[0]=testv[0]; val[1]=testv[1]; val[2]=testv[2]; val[3]=testv[3];
          ffound=0;
          break;
        } else {
          minw=testw[0];
          for (int k=1; k<4; k++) if (testw[k]<minw) minw=testw[k];
          if (minw > minw_best) {
            tbackward=cur_bin->date[i];
            weight[0]=testw[0]; weight[1]=testw[1]; weight[2]=testw[2]; weight[3]=testw[3];
            val[0]=testv[0]; val[1]=testv[1]; val[2]=testv[2]; val[3]=testv[3];
            minw_best=minw;
            ffound=-1;
          }
          //printf("Negative weights!\n");
	}
      //}
    }
    //date[i].write_string(tstring);
    //fprintf(fs, "%24s %10.3f %10.3f %3d\n", tstring, critdist, dmin, minind);
  }

#ifdef PRINT_FLOC
  printf("b: %s %d %d\n", filelist[indf], inds, minind);
#endif

  nchan=data->nchan;

  //fclose(fs);

  return ffound;

}

float amsu_interp_obj::colloc_forward(float lon0,           //lon-lat coords of int. point
                        float lat0,
                        time_class t0,                  //time of interpolation point
                        float critdist,                 //maximum distance to satellite pixel
                        time_class &tforward,            //time of satellite measurement
			long &scanline,
			long &pixnum,
			char *&filename,
			long &nchan,
                        float *&val)                   //four data values
{

//  long np, nchan;

  FILE *fs, *fs2;

  //bin index:
  long lonind, latind;

  long nrec;		//number of records in the bin
  amsu_index_bin *cur_bin;

  //interpolation point:
  long ind, indf, inds;

  //for looking backwards:
  time_class tnext;

  //amsu data:
  amsu_1c_data *data;

  float minw;
  float minw_best;	//highest minimum weight
  float *testv[4];	//equivalent measurements
  float testw[4];

  //distances from each point in the scan line:
  float d;
  float dmin;		//minimum distance
  long minind;		//where in the scan line...

  //for output:
  char tstring[30];
  char tstring1[30];

  int ffound=-1;   //found forwards, found backwards

  //calculate the bin index:
  binind(lon0, lat0, lonind, latind);

  cur_bin=&bin[latind][lonind];
  nrec=bin[latind][lonind].n;

  if (nrec == 0) {
    fprintf(stderr, "Error: index bin (%f, %f) is empty\n", lon0, lat0);
    return ffound;
  }
  if (t0 > cur_bin->date[nrec-1]) {
    t0.write_string(tstring);
    cur_bin->date[nrec-1].write_string(tstring1);  
    fprintf(stderr, "Error: date (%s) later than any in bin (%s)\n", tstring, tstring1);
    return ffound;
  }  

  ind=bin_search(cur_bin->date, nrec, t0, -1)+1;

  minw_best=MINW_EXTRAP;

  //fs=fopen("testint.txt", "w");
  for (long i=ind; i<nrec; i++) {
    indf=cur_bin->file_index[i];
    inds=cur_bin->scan_index[i];

    data=get_data_forward(indf);

    dmin=sdist(data->data[inds].lon[0], data->data[inds].lat[0], lon0, lat0);
    minind=0;
    for (long j=1; j<data->np; j++) {
      d=sdist(data->data[inds].lon[j], data->data[inds].lat[j], lon0, lat0);
      if (d<dmin) {
        dmin=d;
        minind=j;
      }
    }
    //we may have found our point: get the weights, 
    //if none of them are negative, set the flag and break out of the loop
    //since we don't want to do an extrapolation...
    if (dmin < critdist) {
      //printf("dmin=%f; critdist=%f\n", dmin, critdist);
      scanline=inds;
      pixnum=minind;
      filename=filelist[indf];
      tforward=data->data[inds].date;
      val=data->data[inds].bt[minind];
      ffound=dmin;
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

float amsu_interp_obj::colloc_backward(float lon0,           //lon-lat coords of int. point
                        float lat0,
                        time_class t0,                  //time of interpolation point
                        float critdist,                 //maximum distance to satellite pixel
                        time_class &tbackward,            //time of satellite measurement
			long &scanline,
			long &pixnum,
			char *& filename,
			long &nchan,
                        float *& val)                   //four data values
{

//  long np, nchan;

  FILE *fs, *fs2;

  //bin index:
  long lonind, latind;

  long nrec;		//number of records in the bin
  amsu_index_bin *cur_bin;

  //interpolation point:
  long ind, indf, inds;

  //for looking backwards:
  time_class tnext;

  //amsu data:
  amsu_1c_data *data;

  float minw;
  float minw_best;	//highest minimum weight

  //distances from each point in the scan line:
  float d;
  float dmin;		//minimum distance
  long minind;		//where in the scan line...

  //for output:
  char tstring[30];
  char tstring1[30];

  int ffound=-1;   //found forwards, found backwards

  //calculate the bin index:
  binind(lon0, lat0, lonind, latind);

  cur_bin=&bin[latind][lonind];
  nrec=bin[latind][lonind].n;

  if (nrec == 0) {
    fprintf(stderr, "Error: index bin (%f, %f) is empty\n", lon0, lat0);
    return ffound;
  }
  if (t0 > cur_bin->date[nrec-1]) {
    t0.write_string(tstring);
    cur_bin->date[nrec-1].write_string(tstring1);  
    fprintf(stderr, "Error: date (%s) later than any in bin (%s)\n", tstring, tstring1);
    return ffound;
  }  

  ind=bin_search(cur_bin->date, nrec, t0, -1)+1;

  minw_best=MINW_EXTRAP;

  //fs=fopen("testint.txt", "w");
  for (long i=ind; i>=0; i--) {
    indf=cur_bin->file_index[i];
    inds=cur_bin->scan_index[i];

    data=get_data_backward(indf);

    dmin=sdist(data->data[inds].lon[0], data->data[inds].lat[0], lon0, lat0);
    minind=0;
    for (long j=1; j<data->np; j++) {
      d=sdist(data->data[inds].lon[j], data->data[inds].lat[j], lon0, lat0);
      if (d<dmin) {
        dmin=d;
        minind=j;
      }
    }
    //we may have found our point: get the weights, 
    //if none of them are negative, set the flag and break out of the loop
    //since we don't want to do an extrapolation...
    if (dmin < critdist) {
      //printf("dmin=%f; critdist=%f\n", dmin, critdist);
      scanline=inds;
      pixnum=minind;
      filename=filelist[indf];
      tbackward=data->data[inds].date;
      val=data->data[inds].bt[minind];
      ffound=dmin;
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

