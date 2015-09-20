#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "amsu_tinterp_obj.h"
#include "tcoord_defs.h"
#include "global_metric.h"

#define MAX_LL 200
#define MAX_NLOAD 5
#define SWAP_ENDIAN 1

amsu_tinterp_obj::amsu_tinterp_obj(char *indexfile, 		//file containing indexed points
				long mnl, 			//number of data files to keep in memory
				char *basepath, 		//path to data files
				char *nfile, 			//file containing N. hemisphere velocity fields
				char *sfile, 			//file containing S. hemisphere velocity fields
				interpol_index dt1,		//time step
				long nt1):			//number of time steps
		amsu_interp_obj(indexfile, mnl, basepath) {

  FILE *fs;
  char line[MAX_LL];
  long fp_amsustart;		//start of list of amsu files
  long ll;			//line length
  char ** flist1;		//for sorting file-list
  time_class *ts;
  long *sind;			//sorting indices

  //throw-away:
  long np, nchan;

  tion=new traj_int_obj(nfile);
  tios=new traj_int_obj(sfile);

  dt=dt1;
  nt=nt1;

}

amsu_tinterp_obj::~amsu_tinterp_obj() {
  delete tion;
  delete tios;

}

#define MINW_EXTRAP -1	

int amsu_tinterp_obj::int_forward(float lon0,           //lon-lat coords of int. point
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

  float coord0[2], coord1[2], coord2[2];	//trajectory start and end points
  interpol_index tind;
  traj_int_obj *tio;

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

  //integrate the trajectory:
  coord0[0]=lon0;
  coord0[1]=lat0;
  if (lat0 < 0) {
    tio=tios;
  } else {
    tio=tion;
  }
  tind=tio->get_tind(t0);
  tio->integrate(coord0, coord2, tind, dt, nt);

  ind=bin_search(cur_bin->date, nrec, t0, -1)+1;

  minw_best=MINW_EXTRAP;

  //fs=fopen("testint.txt", "w");
  for (long i=ind; i<nrec; i++) {
    indf=cur_bin->file_index[i];
    inds=cur_bin->scan_index[i];

    data=get_data_forward(indf);

    if (tio->get_result(cur_bin->date[i], coord1)!=0) {
      tind+=dt*nt;
      tio->integrate(coord2, coord2, tind, dt, nt);
      tio->get_result(cur_bin->date[i], coord1);
    }

    dmin=sdist(data->data[inds].lon[0], data->data[inds].lat[0], coord1[0], coord1[1]);
    minind=0;
    for (long j=1; j<data->np; j++) {
      d=sdist(data->data[inds].lon[j], data->data[inds].lat[j], coord1[0], coord1[1]);
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
      derive_amsu_int_weights(data, inds, minind, coord1[0], coord1[1], testv, testw);
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

int amsu_tinterp_obj::int_backward(float lon0,           //lon-lat coords of int. point
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

  float coord0[2], coord1[2], coord2[2];	//trajectory start and end points
  interpol_index tind;
  traj_int_obj *tio;

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

  //integrate the trajectory:
  coord0[0]=lon0;
  coord0[1]=lat0;
  if (lat0 < 0) {
    tio=tios;
  } else {
    tio=tion;
  }
  tind=tio->get_tind(t0);
  tio->integrate(coord0, coord2, tind, -dt, nt);

  ind=bin_search(cur_bin->date, nrec, t0, -1)+1;

  //fs=fopen("testint.txt", "w");
  for (long i=ind; i>=0; i--) {
    indf=cur_bin->file_index[i];
    inds=cur_bin->scan_index[i];

    data=get_data_backward(indf);

    if (tio->get_result(cur_bin->date[i], coord1)!=0) {
      tind-=dt*nt;
      tio->integrate(coord2, coord2, tind, -dt, nt);
      tio->get_result(cur_bin->date[i], coord1);
    }

    dmin=sdist(data->data[inds].lon[0], data->data[inds].lat[0], coord1[0], coord1[1]);
    minind=0;
    for (long j=1; j<data->np; j++) {
      d=sdist(data->data[inds].lon[j], data->data[inds].lat[j], coord1[0], coord1[1]);
      if (d<dmin) {
        dmin=d;
        minind=j;
      }
    }

    //printf("%d: dmin=%f\n", i, dmin);
    //we found our point: get the weights, set the flag and break out of the loop:
    if (dmin < critdist) {
      //printf("dmin=%f; critdist=%f\n", dmin, critdist);
      derive_amsu_int_weights(data, inds, minind, coord1[0], coord1[1], testv, testw);
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

float amsu_tinterp_obj::colloc_forward(float lon0,           //lon-lat coords of int. point
                        float lat0,
                        time_class t0,                  //time of interpolation point
                        float critdist,                 //maximum distance to satellite pixel
                        time_class &tforward,            //time of satellite measurement
       		        long &scanline,          //scan number
	                long &pixnum,            //pixel number
        	        char *&filename,         //(pointer to) filename
                	long &nchan,            //total number of channels
	                float *&val)            //(pointer to) values of nearest neighbour
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

  float ffound=-1;   //found forwards, found backwards

  float coord0[2], coord1[2], coord2[2];	//trajectory start and end points
  short hemi;
  interpol_index tind;
  traj_int_obj *tio;

  float lon, lat;

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

  //integrate the trajectory:
  if (lat0 < 0) {
    hemi=-1;
    tio=tios;
  } else {
    hemi=1;
    tio=tion;
  }
  //convert to pole-centred azimuthal equidistant coordinates system:
  tcoord2_2lonlat(coord0[0], coord0[1], -1, hemi, lon0, lat0);
  tind=tio->get_tind(t0);
  tio->integrate(coord0, coord2, tind, dt, nt);

  ind=bin_search(cur_bin->date, nrec, t0, -1)+1;

  minw_best=MINW_EXTRAP;

  //fs=fopen("testint.txt", "w");
  for (long i=ind; i<nrec; i++) {
    indf=cur_bin->file_index[i];
    inds=cur_bin->scan_index[i];

    data=get_data_forward(indf);

    if (tio->get_result(cur_bin->date[i], coord1)!=0) {
      tind+=dt*nt;
      tio->integrate(coord2, coord2, tind, dt, nt);
      tio->get_result(cur_bin->date[i], coord1);
    }
    //convert to regular lon-lat coords:
    tcoord2_2lonlat(coord1[0], coord1[1], 1, hemi, lon, lat);
    //printf("(%f, %f)\n", lon, lat);

    dmin=sdist(data->data[inds].lon[0], data->data[inds].lat[0], lon, lat);
    minind=0;
    for (long j=1; j<data->np; j++) {
      d=sdist(data->data[inds].lon[j], data->data[inds].lat[j], lon, lat);
      if (d<dmin) {
        dmin=d;
        minind=j;
      }
    }

    if (dmin < critdist) {
      //printf("dmin=%f; critdist=%f\n", dmin, critdist);
      tforward=data->data[inds].date;
      scanline=inds;
      pixnum=minind;
      filename=filelist[indf];
      val=data->data[inds].bt[minind];
      ffound=dmin;
      break;
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

float amsu_tinterp_obj::colloc_backward(float lon0,           //lon-lat coords of int. point
                        float lat0,
                        time_class t0,                  //time of interpolation point
                        float critdist,                 //maximum distance to satellite pixel
                        time_class &tbackward,            //time of satellite measurement
       		        long &scanline,          //scan number
	                long &pixnum,            //pixel number
        	        char *&filename,         //(pointer to) filename
                	long &nchan,            //total number of channels
	                float *&val)            //(pointer to) values of nearest neighbour
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

  float ffound=-1;   //found forwards, found backwards

  float coord0[2], coord1[2], coord2[2];	//trajectory start and end points
  short hemi;
  interpol_index tind;
  traj_int_obj *tio;

  float lon, lat;

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

  //integrate the trajectory:
  coord0[0]=lon0;
  coord0[1]=lat0;
  if (lat0 < 0) {
    hemi=-1;
    tio=tios;
  } else {
    hemi=1;
    tio=tion;
  }
  //convert to pole-centred azimuthal equidistant coordinates system:
  tcoord2_2lonlat(coord0[0], coord0[1], -1, hemi, lon0, lat0);
  //printf("(%f, %f); ", lon0, lat0);
  //printf("(%f, %f)\n", coord0[0], coord0[1]);
  tind=tio->get_tind(t0);
  tio->integrate(coord0, coord2, tind, -dt, nt);
  //tio->print_result(stdout);

  ind=bin_search(cur_bin->date, nrec, t0, -1);

  minw_best=MINW_EXTRAP;

  //fs=fopen("testint.txt", "w");
  for (long i=ind; i>=0; i--) {
    indf=cur_bin->file_index[i];
    inds=cur_bin->scan_index[i];

    data=get_data_backward(indf);

    if (tio->get_result(cur_bin->date[i], coord1)!=0) {
      tind-=dt*nt;
      tio->integrate(coord2, coord2, tind, -dt, nt);
      tio->get_result(cur_bin->date[i], coord1);
      //tio->print_result(stdout);
    }
    tcoord2_2lonlat(coord1[0], coord1[1], 1, hemi, lon, lat);
    cur_bin->date[i].write_string(tstring);
    //printf("%s: (%f, %f); ", tstring, coord1[0], coord1[1]);
    //printf("(%f, %f)\n", lon, lat);

    dmin=sdist(data->data[inds].lon[0], data->data[inds].lat[0], lon, lat);
    minind=0;
    for (long j=1; j<data->np; j++) {
      d=sdist(data->data[inds].lon[j], data->data[inds].lat[j], lon, lat);
      if (d<dmin) {
        dmin=d;
        minind=j;
      }
    }

    if (dmin < critdist) {
      //printf("dmin=%f; critdist=%f\n", dmin, critdist);
      tbackward=data->data[inds].date;
      scanline=inds;
      pixnum=minind;
      filename=filelist[indf];
      val=data->data[inds].bt[minind];
      ffound=dmin;
      break;
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

