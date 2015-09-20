#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "peteys_tmpl_lib.h"
#include "global_metric.h"
#include "read_amsu.h"

#define DTOL 0.001

amsu_1c_data * interpolate_amsu_1c_scan (amsu_1c_data * data1, long npnew, long offset) {
  amsu_1c_data *result;
  float lon1, lon2;
  float lat1, lat2;
  float *bt1, *bt2;

  double *intcoef;		//interpolation coefficient
  long ind;
  double frac;

  result=new amsu_1c_data;
  result->data=new amsu_1c_rec[data1->nscan];
  result->nscan=data1->nscan;
  result->np=npnew;
  result->nchan=data1->nchan;
  result->missing=data1->missing;

  intcoef=new double[npnew];
  for (long i=0; i<npnew; i++) {
    intcoef[i]=(i+offset)*((double) data1->np-1)/((double) npnew+2.*offset-1);
    //no extrapolation:
    if (intcoef[i] >= data1->np-1) intcoef[i]=data1->np-1;
    		else if (intcoef[i] < 0) intcoef[i]=0;
    //printf("%f ", intcoef[i]);
  }
  //printf("\n");

  //start the interpolation:
  for (long i=0; i<data1->nscan; i++) {
    //allocate the new data:
    result->data[i].lon=new float[result->np];
    result->data[i].lat=new float[result->np];
    result->data[i].bt=new float *[result->np];

    //first, copy the date:
    result->data[i].date=data1->data[i].date;

    for (long j=0; j<result->np; j++) {
      ind=(long) intcoef[j];
      if (ind >= data1->np-1) ind=data1->np-2;
      else if (ind < 0) ind=0;
      frac=intcoef[j]-(double) ind;

      //interpolate coordinates:      
      interpolate_lon_lat(data1->data[i].lon[ind], data1->data[i].lat[ind], 
		      		data1->data[i].lon[ind+1], data1->data[i].lat[ind+1], frac, 
				result->data[i].lon[j], result->data[i].lat[j]);

      //interpolate brightness temperatures:
      //
      //no extrapolation:
      if (ind >= data1->np-1) {
        ind=data1->np-2;
	frac=1;
      } else if (ind < 0) {
        ind=0;
	frac=0;
      }
      
      result->data[i].bt[j]=new float[data1->nchan];
      for (long k=0; k<result->nchan; k++) {
        if (data1->data[i].bt[ind][k] <= data1->missing ||
			data1->data[i].bt[ind+1][k] <= data1->missing) {
          result->data[i].bt[j][k]=data1->missing;
	} else {
          result->data[i].bt[j][k]=(1-frac)*data1->data[i].bt[ind][k]+
			frac*data1->data[i].bt[ind+1][k];
	}
      }
    }

  }

  delete [] intcoef;

  return result;

}


//interpolates the points on the scan line in data1 to the equivalent of data2;
//data2 must have the same time grids as data1...
amsu_1c_data * interpolate_amsu_1c_scan2 (amsu_1c_data * data1, amsu_1c_data *data2) {
  amsu_1c_data *result;
  float lon1, lon2;
  float lat1, lat2;
  float *bt1, *bt2;

  float *s1, *s2;		//distance along the scan line
  float d;

  double intcoef;		//interpolation coefficient
  double frac;			//fractional components
  long ind;			//interpolation index

  assert(data1->nscan == data2->nscan);

  result=new amsu_1c_data;
  result->data=new amsu_1c_rec[data1->nscan];
  result->nscan=data1->nscan;
  result->np=data2->np;
  result->nchan=data1->nchan;
  result->missing=data1->missing;

  s1=new float[data1->np];
  s2=new float[data2->np];

  //start the interpolation:
  for (long i=0; i<data1->nscan; i++) {
    //allocate the new data:
    result->data[i].lon=new float[result->np];
    result->data[i].lat=new float[result->np];
    result->data[i].bt=new float *[result->np];

    //first, copy the date:
    result->data[i].date=data1->data[i].date;

    //calculate distances along the scan line:
    s1[0]=0;
    for (long j=1; j<data1->np; j++) s1[j]=s1[j-1]+
		    sqrt(sdist(data1->data[i].lon[j-1], data1->data[i].lat[j-1], 
		    	data1->data[i].lon[j], data1->data[i].lat[j]));
    for (long j=0; j<data1->np; j++) printf("%8f ", s1[j]);
    printf("\n");

    //find offset of one scan line from the other:
    s2[0]=sqrt(sdist(data1->data[i].lon[0], data1->data[i].lat[0], 
			    data2->data[i].lon[0], data2->data[i].lat[0]));
    
    //is the offset positive or negative?:
    d=sqrt(sdist(data1->data[i].lon[1], data1->data[i].lat[1], 
			    data2->data[i].lon[0], data2->data[i].lat[0]));
    if (d+s2[0]-DTOL > s1[1]-s1[0]) s2[0]=-s2[0];

    for (long j=1; j<data2->np; j++) s2[j]=s2[j-1]+
		    sqrt(sdist(data2->data[i].lon[j-1], data2->data[i].lat[j-1], 
		    	data2->data[i].lon[j], data2->data[i].lat[j]));
    for (long j=0; j<data2->np; j++) printf("%8f ", s2[j]);
    printf("\n");

    for (long j=0; j<data2->np; j++) {
      //simply copy longitudes and latitudes:
      result->data[i].lon[j]=data2->data[i].lon[j];
      result->data[i].lat[j]=data2->data[i].lat[j];

      //calculate interpolation coefficients:      
      intcoef=interpolate(s1, data1->np, s2[j], -1);
      printf("%8.2f ", intcoef);
      ind=(long) intcoef;
      if (ind >= data1->np-1) ind=data1->np-2;
      else if (ind < 0) ind=0;
      frac=intcoef-(double) ind;

      //interpolate brightness temperatures:
      result->data[i].bt[j]=new float[data1->nchan];
      for (long k=0; k<result->nchan; k++) {
        if (data1->data[i].bt[ind][k] <= data1->missing ||
			data1->data[i].bt[ind+1][k] <= data1->missing) {
          result->data[i].bt[j][k]=data1->missing;
	} else {
          result->data[i].bt[j][k]=(1-frac)*data1->data[i].bt[ind][k]+
			frac*data1->data[i].bt[ind+1][k];
	}
      }
    }
    printf("\n");

  }

  return result;

}

//copies an amsu_1c_data structure, but only the date part:
amsu_1c_data * cp_amsu_1c_date(amsu_1c_data * data1) {
  amsu_1c_data *data2;

  data2=new amsu_1c_data;  
  data2->nscan=data1->nscan;
  data2->np=0;
  data2->nchan=0;
  data2->missing=data1->missing;

  data2->data=new amsu_1c_rec[data2->nscan];

  for (long i=0; i<data2->nscan; i++) {
    data2->data[i].date=data1->data[i].date;
    data2->data[i].lon=NULL;
    data2->data[i].lat=NULL;
    data2->data[i].bt=NULL;
  }

  return data2;

}

//interpolates amsu 1c data in time:
//*** note: data2 has all the scan-lines and dates for each scan-line,
//but is otherwise empty...

void interpolate_amsu_1c_date(amsu_1c_data * data1, amsu_1c_data *data2) {
  time_class date[data1->nscan];
  double ind;
  long lind=-1;
  double frac;
  long last_ind=-1;
  float *lon1, *lon2;
  float *lat1, *lat2;
  float **bt1, **bt2;

  data2->np=data1->np;
  data2->nchan=data1->nchan;
  data2->missing=data1->missing;

  for (long i=0; i<data1->nscan; i++) date[i]=data1->data[i].date;
  for (long i=0; i<data2->nscan; i++) {
    ind=interpolate(date, data1->nscan, data2->data[i].date, last_ind);
    lind=(long) ind;
    last_ind=lind;
    if (lind < 0) lind=0; else if (lind >= data1->nscan-1) lind=data1->nscan-2;

    //set the temporaries for each side of the grid interval:
    lon1=data1->data[lind].lon;
    lon2=data1->data[lind+1].lon;
    lat1=data1->data[lind].lat;
    lat2=data1->data[lind+1].lat;
    bt1=data1->data[lind].bt;
    bt2=data1->data[lind+1].bt;

    //allocate the variables:
    data2->data[i].lon=new float[data2->np];
    data2->data[i].lat=new float[data2->np];
    data2->data[i].bt=new float*[data2->np];

    //do the interpolation:
    frac=ind-(double) lind;
    for (long j=0; j<data2->np; j++) {
      interpolate_lon_lat(lon1[j], lat1[j], lon2[j], lat2[j], frac, 
		      data2->data[i].lon[j], data2->data[i].lat[j]);
      data2->data[i].bt[j]=new float[data2->nchan];
      for (long k=0; k<data2->nchan; k++) {
        if (bt1[j][k] <= data1->missing || bt2[j][k] <= data1->missing) {
          data2->data[i].bt[j][k]=data1->missing;
        } else {
          data2->data[i].bt[j][k]=(1-frac)*bt1[j][k]+frac*bt2[j][k];
        }
      }
    }

  }
}

//reads in limb-correction-coefficients from a text file:
float ***read_amsu_lcc(char *filename, long np, long nchan) {
  FILE *fs;
  float ***c;
  long err;

  fs=fopen(filename, "r");
  if (fs == NULL) return NULL;

  c=new float **[nchan];
  for (long i=0; i<nchan; i++) {
    c[i]=new float *[np];
    c[i][0]=new float[np*(nchan+1)];
    for (long j=0; j<np; j++) {
      c[i][j]=c[i][0]+j*(nchan+1);
      for (long k=0; k<=nchan; k++) {
        err=fscanf(fs, "%f", c[i][j]+k);
        if (err!=1) {
          fprintf(stderr, "read_amsu_lcc: Error reading file, %s\n", filename);
          exit(-10);
        }
      }
    }
  }

  return c;
}

//deletes limb-correction-coefficients:
void delete_amsu_lcc(float ***c, long np) {
  for (long i=0; i<np; i++) {
    delete [] c[i][0];
    delete [] c[i];
  }

  delete [] c;
}

//applies limb-correction-coefficients:
amsu_1c_data * amsu_lc(amsu_1c_data *data, float ***c) {
  int missflag;

  amsu_1c_data *result;

  result=cp_amsu_1c_date(data);

  result->nchan=data->nchan;
  result->np=data->np;
  result->missing=data->missing;
  
  for (long i=0; i<data->nscan; i++) {
    //allocate data records:
    result->data[i].lon=new float[result->np];
    result->data[i].lat=new float[result->np];
    result->data[i].bt=new float*[result->np];

    for (long j=0; j<data->np; j++) {
      //copy positional data:
      result->data[i].lon[j]=data->data[i].lon[j];
      result->data[i].lat[j]=data->data[i].lat[j];
      //allocate channel bts:
      result->data[i].bt[j]=new float[result->nchan];

      //check for missing data:
      missflag=0;
      for (long k=0; k<data->nchan; k++) {
        if (data->data[i].bt[j][k]<=data->missing) {
          missflag=1;
          break;
        }
      }
      if (missflag) {
        //if one channel in the pixel has missing data,
	//assign a missing value to all the channels in the result:
        for (long k=0; k<data->nchan; k++) {
          result->data[i].bt[j][k]=result->missing;
	}
      } else {
        //otherwise, apply the coefficients:
        for (long k=0; k<data->nchan; k++) {
	  //first element has constant term:
	  result->data[i].bt[j][k]=c[k][j][0];
	  for (long m=0; m<data->nchan; m++) {
	    result->data[i].bt[j][k]+=data->data[i].bt[j][m]*c[k][j][m+1];
	  }
	}
      }
    }
  }

  return result;
}
      
