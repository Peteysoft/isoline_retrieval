#include <stdio.h>
#include <math.h>

#include "three_d_fields.h"

float *** allocate_3d_field(long nlev, long nlat, long nlon) {
  FILE *fs;
  float *** data;

  data=new float **[nlev];
  data[0]=new float *[nlat];
  data[0][0]=new float [nlev*nlat*nlon];
  for (long j=1; j<nlat; j++) {
    data[0][j]=&data[0][0][nlon*j];
  }
  for (long i=1; i<nlev; i++) {
    data[i]=new float *[nlat];
    for (long j=0; j<nlat; j++) {
      data[i][j]=&data[0][0][nlon*(j+nlat*i)];
    }
  }

  return data;

}

void delete_3d_field(float *** field, long nlev, long nlat, long nlon) {
  for (long i=1; i<nlev; i++) delete [] field[i];
  delete [] field [0][0];
  delete [] field [0];
  delete [] field;
}

float *** read_grads_bin(char *filename, long nlev, long nlat, long nlon) {
  FILE *fs;
  float *** data;
  long nread;

  fs=fopen(filename, "r");
  if (fs == NULL) {
    printf("File %s not found\n", filename);
    return NULL;
  }

  data=allocate_3d_field(nlev, nlat, nlon);

  nread=fread(data[0][0], sizeof(float), nlev*nlat*nlon, fs);
  if (nread != nlev*nlat*nlon) {
    printf("File read error: %s\n", filename);
    delete_3d_field(data, nlev, nlat, nlon);
    data = NULL;
  }

  fclose(fs);

  return data;

}

int read_grads_bin(char *filename, long nlev, long nlat, long nlon, 
		float *** data) {
  FILE *fs;
  long nread;

  fs=fopen(filename, "r");
  if (fs == NULL) {
    printf("File %s not found\n", filename);
    return -1;
  }

  nread=fread(data[0][0], sizeof(float), nlev*nlat*nlon, fs);
  if (nread != nlev*nlat*nlon) {
    printf("File read error: %s\n", filename);
  }

  fclose(fs);

  return nread;

}

float *** calc_height(float *** pp, float *** tt, float *** vmr, 
		long nlev, long nlat, long nlon) {

  float *** zz;
  float wvp;		//water vapour pressure
  float Req1, Req2;	//equivalent gas constant
  float dz;		//layer thickness
  float c;		//constant of integration

  zz=allocate_3d_field(nlev, nlat, nlon);

  for (long i=0; i<nlon; i++) {
    for (long j=0; j<nlat; j++) {
      zz[nlev-1][j][i]=0;		//altitude at bottom is 0, right?
      wvp=vmr[nlev-1][j][i]*pp[nlev-1][j][i];	//calculate water vapour pressure from vmr

      //calculate the equivalent gas constant:
      Req1=pp[nlev-1][j][i]/((pp[nlev-1][j][i]-wvp)/RDRYAIR+wvp/RWATERVAPOUR);
      for (long k=nlev-2; k>=0; k--) {
        wvp=vmr[k][j][i]*pp[k][j][i];
        Req2=pp[k][j][i]/((pp[k][j][i]-wvp)/RDRYAIR+wvp/RWATERVAPOUR);
        //constant of integration:
        c=(Req1*tt[k+1][j][i]+Req2*tt[k][j][i])/GACC/2;
        //calculate layer thickness:
        dz=c*log(pp[k+1][j][i]/pp[k][j][i]);
        zz[k][j][i]=zz[k+1][j][i]+dz;
        Req1=Req2;
      }
    }
  }

  return zz;

}


