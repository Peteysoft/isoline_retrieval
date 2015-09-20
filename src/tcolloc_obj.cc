#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tcolloc_obj.h"
#include "tcoord_defs.h"
#include "global_metric.h"

#define MAX_LL 200
#define MAX_NLOAD 5
#define SWAP_ENDIAN 1

tcolloc_obj::tcolloc_obj(char *indexfile, 		//file containing indexed points
				long mnl, 			//number of data files to keep in memory
				char *basepath, 		//path to data files
				char *nfile, 			//file containing N. hemisphere velocity fields
				char *sfile, 			//file containing S. hemisphere velocity fields
				interpol_index dt1,		//time step
				long nt1):			//number of time steps
		colloc_obj(indexfile, mnl, basepath) {

  tion=new traj_int_obj(nfile);
  tios=new traj_int_obj(sfile);

  dt=dt1;
  nt=nt1;

}

tcolloc_obj::~tcolloc_obj() {
  delete tion;
  delete tios;

}

//set the page size for the velocity field reading routines:
void tcolloc_obj::set_page_size(long page_size) {
  tion->set_page_size(page_size);
  tios->set_page_size(page_size);
}


float tcolloc_obj::colloc_forward(float lon0,           //lon-lat coords of int. point
                        float lat0,
                        time_class t0,                  //time of interpolation point
                        float critdist,                 //maximum distance to satellite pixel
                        time_class &tforward,            //time of satellite measurement
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
  gome_data *data;

  float minw;
  float minw_best;	//highest minimum weight

  //distances from each point in the scan line:
  float d;
  float dmin;		//minimum distance
  long minind;		//where in the scan line...
  int4_t minsind;
  long minfind;

  //for output:
  char tstring[30];
  char tstring1[30];

  float ffound=-1;   //found forwards, found backwards

  float coord0[2], coord1[2], coord2[2];	//trajectory start and end points
  short hemi;
  interpol_index tind;
  traj_int_obj *tio;

  float lon, lat;

  double tdiff;

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

  t0.write_string(tstring1);
  //printf("(%f, %f); %s\n", lon0, lat0, tstring);

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

  //we want to check our results:
  //if (hemi > 0) tio->print_result_N(stdout); else tio->print_result_S(stdout);

  ind=bin_search(cur_bin->date, nrec, t0, -1)+1;
  if (ind >= nrec) return ffound;

  indf=cur_bin->file_index[ind];
  inds=cur_bin->scan_index[ind];

  data=get_data_forward(indf);

  tio->get_result(cur_bin->date[ind], coord1);
  //convert to regular lon-lat coords:
  tcoord2_2lonlat(coord1[0], coord1[1], 1, hemi, lon, lat);

  tdiff=(tind-tio->get_tind(cur_bin->date[ind]))/dt/nt;

  dmin=sdist(data->data[inds].lon, data->data[inds].lat, lon, lat)+tdiff*tdiff*critdist;
  minfind=indf;
  minsind=inds;

  //fs=fopen("testint.txt", "w");
  for (long i=ind+1; i<nrec; i++) {
    indf=cur_bin->file_index[i];
    inds=cur_bin->scan_index[i];

    data=get_data_forward(indf);

    if (tio->get_result(cur_bin->date[i], coord1)!=0) {
      break;
      //tind+=dt*nt;
      //tio->integrate(coord2, coord2, tind, dt, nt);
      //tio->get_result(cur_bin->date[i], coord1);
    }
    //convert to regular lon-lat coords:
    tcoord2_2lonlat(coord1[0], coord1[1], 1, hemi, lon, lat);

    tdiff=(tind-tio->get_tind(cur_bin->date[i]))/dt/nt;

    d=sdist(data->data[inds].lon, data->data[inds].lat, lon, lat)+tdiff*tdiff*critdist;
    cur_bin->date[i].write_string(tstring);
    //printf("%s; |(%f, %f), (%f, %f)|=%f\n", tstring, data->data[inds].lon, data->data[inds].lat, lon, lat, d);

    //printf("%f, %s, %s\n", sqrt(tdiff*tdiff*critdist), tstring1, tstring);

    //printf("(%f, %f) - (%f, %f); d=%f\n", lon, lat, 
	//	data->data[inds].lon, data->data[inds].lat, sqrt(d));

    if (d<dmin) {
      dmin=d;
      minfind=indf;
      minsind=inds;
      traj_lon=lon;
      traj_lat=lat;
    }

    //fprintf(fs, "%24s %10.3f %10.3f %3d\n", tstring, critdist, dmin, minind);
    //printf("%24s %10.3f %10.3f %3d\n", tstring, critdist, dmin, minind);
  }

#ifdef PRINT_FLOC
  printf("f: %s %d %d\n", filelist[indf], inds, minind);
#endif

  data=get_data_backward(minfind);

  //printf("dmin=%f; critdist=%f\n", dmin, critdist);
  tforward=data->data[minsind].date;
  pixnum=minsind;
  filename=filelist[minfind];
  val=data->data[minsind].counts;

  tdiff=(tind-tio->get_tind(tforward))/dt/nt;
  //printf("%f %f %f %f\n", tdiff, tind-tio->get_tind(tforward), critdist, dt*nt);
  tforward.write_string(tstring);
  //printf("%s: (%f, %f)\n", tstring, data->data[minsind].lon, data->data[minsind].lat);
  ffound=dmin-tdiff*tdiff*critdist;

  nchan=data->nchan;

  //fclose(fs);

  return ffound;

}

float tcolloc_obj::colloc_backward(float lon0,           //lon-lat coords of int. point
                        float lat0,
                        time_class t0,                  //time of interpolation point
                        float critdist,                 //maximum distance to satellite pixel
                        time_class &tbackward,            //time of satellite measurement
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
  gome_data *data;

  float minw;
  float minw_best;	//highest minimum weight
  int4_t minsind;
  long minfind;

  //distances from each point in the scan line:
  float d;
  float dmin;		//minimum distance

  //for output:
  char tstring[30];
  char tstring1[30];

  float ffound=-1;   //found forwards, found backwards

  float coord0[2], coord1[2], coord2[2];	//trajectory start and end points
  short hemi;
  interpol_index tind;
  traj_int_obj *tio;

  float lon, lat;

  double tdiff;

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
  if (ind == -1) return -1;

  indf=cur_bin->file_index[ind];
  inds=cur_bin->scan_index[ind];

  data=get_data_backward(indf);

  tio->get_result(cur_bin->date[ind], coord1);
  tcoord2_2lonlat(coord1[0], coord1[1], 1, hemi, lon, lat);

  tdiff=(tind-tio->get_tind(cur_bin->date[ind]))/dt/nt;
  dmin=sdist(data->data[inds].lon, data->data[inds].lat, lon, lat)+tdiff*tdiff*critdist;
  minfind=indf;
  minsind=inds;

  //fs=fopen("testint.txt", "w");
  for (long i=ind-1; i>=0; i--) {
    indf=cur_bin->file_index[i];
    inds=cur_bin->scan_index[i];

    data=get_data_backward(indf);

    if (tio->get_result(cur_bin->date[i], coord1)!=0) {
      break;
    //  tind-=dt*nt;
    //  tio->integrate(coord2, coord2, tind, -dt, nt);
    //  tio->get_result(cur_bin->date[i], coord1);
      //tio->print_result(stdout);
    }

    tcoord2_2lonlat(coord1[0], coord1[1], 1, hemi, lon, lat);
    cur_bin->date[i].write_string(tstring);
    //printf("%s: (%f, %f); ", tstring, coord1[0], coord1[1]);
    //printf("(%f, %f)\n", lon, lat);

    tdiff=(tind-tio->get_tind(cur_bin->date[i]))/dt/nt;
    d=sdist(data->data[inds].lon, data->data[inds].lat, lon, lat)+tdiff*tdiff*critdist;
    minfind=indf;
    minsind=inds;

    if (d < dmin) {
      dmin=d;
      minfind=indf;
      minsind=inds;
      traj_lon=lon;
      traj_lat=lat;
    }

    //date[i].write_string(tstring);
    //fprintf(fs, "%24s %10.3f %10.3f %3d\n", tstring, critdist, dmin, minind);
  }

#ifdef PRINT_FLOC
  printf("f: %s %d %d\n", filelist[indf], inds, minind);
#endif

  data=get_data_forward(minfind);

  //printf("dmin=%f; critdist=%f\n", dmin, critdist);
  tbackward=data->data[minsind].date;
  pixnum=minsind;
  filename=filelist[minfind];
  val=data->data[minsind].counts;

  tdiff=(tind-tio->get_tind(tbackward))/dt/nt;
  //printf("%f\n", tdiff);
  ffound=dmin-tdiff*tdiff*critdist;

  nchan=data->nchan;

  //fclose(fs);

  return ffound;

}

