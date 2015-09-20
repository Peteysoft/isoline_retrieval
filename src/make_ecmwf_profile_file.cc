#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "time_class.h"

#include "three_d_fields.h"
#include "ecmwf_sampling2.h"

#define NLON 240
#define NLAT 121

void get_ecmwf_field(time_class date, char *var, float ***&data) {
  char filename[200];
  short yy, mm, dd, hh, min;
  float s;

  date.get_fields(yy, mm, dd, hh, min, s);

  sprintf(filename, "%s/%4.4d/%2.2d/%2.2d/%2.2d/%2s%4.4d%2.2d%2.2d%2.2d", 
		  "/equinox/ecmwf", yy, mm, dd, hh, 
		var, yy, mm, dd, hh);

  printf("Reading file: %s\n", filename);

  if (data == NULL) data=read_grads_bin(filename, NLEVELS, NLAT, NLON);
  		else read_grads_bin(filename, NLEVELS, NLAT, NLON, data);

}


int main (int argc, char **argv) {
  FILE *fs;
  FILE *outfs;
  char tstring[20];
  time_class t, told;
  long lonind, latind;
  int ch;
  long n, k;

  float ***tt=NULL;
  float ***pp=NULL;
  float ***q=NULL;
  float ***cl=NULL;
  float ***ci=NULL;
  float ***cc=NULL;
  float ***uu=NULL;
  float ***vv=NULL;

  float req1, req2;	//equivalent gas constant
  float wvp;		//water vapour partial-pressure
  float c;		//constant of integration
  float z1;		//height at first level

  float frac;

  ecmwf_profile prof;

  if (argc != 3) {
    printf("syntax: make_ecmwf_profile_file profile_list outfile\n");
    printf("	where:\n");
    printf("profile_list:	file with list of dates and coordinates\n");
    printf("outfile:		binary output file\n");
    exit(1);
  }

  fs=fopen(argv[1], "r");

  //count the number of lines in the file:
  n=0;  k=0;  c=0;
  while (ch != EOF) {
    ch=fgetc(fs);
    if (ch=='\n') {
      if (k!=0) n++;
      k=0;
    }
    k++;
  }

  outfs=fopen(argv[2], "w");

  //read in all the data:
  rewind(fs);
  fscanf(fs, "%15s%5d%5d", tstring, &lonind, &latind);
  told.read_string(tstring);

  get_ecmwf_field(told, "PP", pp);
  get_ecmwf_field(told, "TT", tt);
  get_ecmwf_field(told, "SQ", q);

  get_ecmwf_field(told, "CL", cl);
  get_ecmwf_field(told, "CI", ci);
  get_ecmwf_field(told, "CC", cc);

  get_ecmwf_field(told, "UU", uu);
  get_ecmwf_field(told, "VV", vv);

  //calculate the thickness of the first layer:
  wvp=q[NLEVELS-1][latind][lonind]*pp[NLEVELS-1][latind][lonind]*1.60627;
  req1=pp[NLEVELS-1][latind][lonind]/((pp[NLEVELS-1][latind][lonind]-wvp)/RDRYAIR+
		wvp/RWATERVAPOUR);
  wvp=q[NLEVELS-2][latind][lonind]*pp[NLEVELS-2][latind][lonind]*1.60627;
  req2=pp[NLEVELS-2][latind][lonind]/((pp[NLEVELS-2][latind][lonind]-wvp)/RDRYAIR+
		wvp/RWATERVAPOUR);
  c=(req1*tt[NLEVELS-1][latind][lonind]+req2*tt[NLEVELS-2][latind][lonind])/GACC/2;
  z1=c*log(pp[NLEVELS-1][latind][lonind]/pp[NLEVELS-2][latind][lonind]); 

  //insert location:
  prof.date=told;
  prof.lon=lonind*360./NLON;
  prof.lat=latind*180./(NLAT-1)-90;

  //insert surface variables:
  prof.p0=pp[NLEVELS-1][latind][lonind];
  prof.p0*=100;
  prof.t0=tt[NLEVELS-1][latind][lonind];

  //calculate two meter variables:
  frac=2./z1;
  prof.p2m=pp[NLEVELS-1][latind][lonind]*(1-frac)+pp[NLEVELS-2][latind][lonind]*frac;
  prof.p2m*=100;
  prof.t2m=tt[NLEVELS-1][latind][lonind]*(1-frac)+tt[NLEVELS-2][latind][lonind]*frac;
  prof.q2m=q[NLEVELS-1][latind][lonind]*(1-frac)+q[NLEVELS-2][latind][lonind]*frac;
//  prof.q2m*=1.60627;

  prof.u2m=uu[NLEVELS-1][latind][lonind]*(1-frac)+uu[NLEVELS-2][latind][lonind]*frac;
  prof.v2m=vv[NLEVELS-1][latind][lonind]*(1-frac)+vv[NLEVELS-2][latind][lonind]*frac;

  //insert the rest of the profile:
  for (long i=0; i<NLEVELS; i++) {
    prof.profile[i].t=tt[i][latind][lonind];
    prof.profile[i].q=q[i][latind][lonind];	//*1.60627;

    prof.profile[i].cl=cl[i][latind][lonind];
    prof.profile[i].ci=ci[i][latind][lonind];
    prof.profile[i].cc=cc[i][latind][lonind];
  }

  write_ecmwf_profile(outfs, &prof);

  for (long i=1; i<n; i++) {
    fscanf(fs, "%15s%5d%5d", tstring, &lonind, &latind);
    t.read_string(tstring);

    if (t!=told) {
      told=t;
      get_ecmwf_field(told, "PP", pp);
      get_ecmwf_field(told, "TT", tt);
      get_ecmwf_field(told, "SQ", q);

      get_ecmwf_field(told, "CL", cl);
      get_ecmwf_field(told, "CI", ci);
      get_ecmwf_field(told, "CC", cc);

      get_ecmwf_field(told, "UU", uu);
      get_ecmwf_field(told, "VV", vv);
    }

    //calculate the thickness of the first layer:
    wvp=q[NLEVELS-1][latind][lonind]*pp[NLEVELS-1][latind][lonind]*1.60627;
    req1=pp[NLEVELS-1][latind][lonind]/((pp[NLEVELS-1][latind][lonind]-wvp)/RDRYAIR+
		wvp/RWATERVAPOUR);
    wvp=q[NLEVELS-2][latind][lonind]*pp[NLEVELS-2][latind][lonind]*1.60627;
    req2=pp[NLEVELS-2][latind][lonind]/((pp[NLEVELS-2][latind][lonind]-wvp)/RDRYAIR+
		wvp/RWATERVAPOUR);
    c=(req1*tt[NLEVELS-1][latind][lonind]+req2*tt[NLEVELS-2][latind][lonind])/GACC/2;
    z1=c*log(pp[NLEVELS-1][latind][lonind]/pp[NLEVELS-2][latind][lonind]); 

    //insert location:
    prof.date=told;
    prof.lon=lonind*360./NLON-180;
    prof.lat=90-latind*180./(NLAT-1);

    //insert surface variables:
    prof.p0=pp[NLEVELS-1][latind][lonind]*100;
    prof.t0=tt[NLEVELS-1][latind][lonind];

    //calculate two meter variables:
    frac=2./z1;
    prof.p2m=pp[NLEVELS-1][latind][lonind]*(1-frac)+pp[NLEVELS-2][latind][lonind]*frac;
    prof.p2m*=100;
    prof.t2m=tt[NLEVELS-1][latind][lonind]*(1-frac)+tt[NLEVELS-2][latind][lonind]*frac;
    prof.q2m=q[NLEVELS-1][latind][lonind]*(1-frac)+q[NLEVELS-2][latind][lonind]*frac;
    //prof.q2m*=1.60627;

    prof.u2m=uu[NLEVELS-1][latind][lonind]*(1-frac)+uu[NLEVELS-2][latind][lonind]*frac;
    prof.v2m=vv[NLEVELS-1][latind][lonind]*(1-frac)+vv[NLEVELS-2][latind][lonind]*frac;

    //insert the rest of the profile:
    for (long j=0; j<NLEVELS; j++) {
      prof.profile[j].t=tt[j][latind][lonind];
      prof.profile[j].q=q[j][latind][lonind];	//*1.60627;

      prof.profile[j].cl=cl[j][latind][lonind];
      prof.profile[j].ci=ci[j][latind][lonind];
      prof.profile[j].cc=cc[j][latind][lonind];
    }

    write_ecmwf_profile(outfs, &prof);
  }

  fclose(fs);

  fclose(outfs);

  delete_3d_field(pp, NLEVELS, NLAT, NLON);
  delete_3d_field(tt, NLEVELS, NLAT, NLON);
  delete_3d_field(q, NLEVELS, NLAT, NLON);

}

