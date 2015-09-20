#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "peteys_tmpl_lib.h"
#include "amsu_interp_base.h"

#include "global_metric.h"

#define MAX_LL 200
#define SWAP_ENDIAN 1

#define COMMENT_LENGTH 200

//#define CHECK_INDEX_FILE 

//#define PRINT_FLOC

//once a point has been located, find three others and compute the weights:
void derive_amsu_int_weights(amsu_1c_data *data,//amsu data
                        long inds,              //the scan index
                        long indp,              //swath index
                        float lon,              //longitude of point
                        float lat,              //latitude of the point
                        float *val[4],           //values at each point
                        float weight[4]) {     //returned weights

  float d;
  float d1, d2, d3, d4;         //distances to each of the four points
  float s;
  float s1, s2, s3, s4;		//side lengths
  long inds1, indp1;           //indices of adjacent points

  d1=sdist(lon, lat, data->data[inds].lon[indp], data->data[inds].lat[indp]);

  //get four pixels that enclose the test point (if possible):
  if (inds == data->nscan-1) {
    inds1=inds-1;
    d2=sdist(lon, lat, data->data[inds1].lon[indp], data->data[inds1].lat[indp]);
    s1=sdist(data->data[inds].lon[indp], data->data[inds].lat[indp],
		  data->data[inds1].lon[indp], data->data[inds1].lat[indp]);
  } else if (inds==0) {
    inds1=1;
    d2=sdist(lon, lat, data->data[inds1].lon[indp], data->data[inds1].lat[indp]);
    s1=sdist(data->data[inds].lon[indp], data->data[inds].lat[indp],
		  data->data[inds1].lon[indp], data->data[inds1].lat[indp]);
  } else {
    //check distance:
    d2=sdist(lon, lat, data->data[inds+1].lon[indp], data->data[inds+1].lat[indp]);
    s1=sdist(data->data[inds].lon[indp], data->data[inds].lat[indp],
		  data->data[inds+1].lon[indp], data->data[inds+1].lat[indp]);
    d=sdist(lon, lat, data->data[inds-1].lon[indp], data->data[inds-1].lat[indp]);
    s=sdist(data->data[inds].lon[indp], data->data[inds].lat[indp],
		  data->data[inds-1].lon[indp], data->data[inds-1].lat[indp]);
    if (d1-d2+s1 >= 0) {
      inds1=inds+1;
    } else {
      inds1=inds-1;
      d2=d;
      s1=s;
    }
  }

  if (indp == data->np-1) {
    indp1=indp-1;
    d4=sdist(lon, lat, data->data[inds].lon[indp1], data->data[inds].lat[indp1]);
    s4=sdist(data->data[inds].lon[indp], data->data[inds].lat[indp],
		  data->data[inds].lon[indp1], data->data[inds].lat[indp1]);
  } else if (indp == 0) {
    indp1=1;
    d4=sdist(lon, lat, data->data[inds].lon[indp1], data->data[inds].lat[indp1]);
    s4=sdist(data->data[inds].lon[indp], data->data[inds].lat[indp],
		  data->data[inds].lon[indp1], data->data[inds].lat[indp1]);
  } else {
    d4=sdist(lon, lat, data->data[inds].lon[indp+1], data->data[inds].lat[indp+1]);
    s4=sdist(data->data[inds].lon[indp], data->data[inds].lat[indp],
		  data->data[inds].lon[indp+1], data->data[inds].lat[indp+1]);
    d=sdist(lon, lat, data->data[inds].lon[indp-1], data->data[inds].lat[indp-1]);
    s=sdist(data->data[inds].lon[indp], data->data[inds].lat[indp],
		  data->data[inds].lon[indp-1], data->data[inds].lat[indp-1]);
    if (d1-d4+s4 >= 0) {
      indp1=indp+1;
    } else {
      indp1=indp-1;
      d4=d;
      s4=s;
    }
  }

  d3=sdist(lon, lat, data->data[inds1].lon[indp1], data->data[inds1].lat[indp1]);
/*
  printf("\n(%10.2f, %10.2f)\n\n", lon, lat);
  printf("(%10.2f, %10.2f)\n", data->data[inds].lon[indp], data->data[inds].lat[indp]);
  printf("(%10.2f, %10.2f)\n", data->data[inds1].lon[indp], data->data[inds1].lat[indp]);
  printf("(%10.2f, %10.2f)\n", data->data[inds1].lon[indp1], data->data[inds1].lat[indp1]);
  printf("(%10.2f, %10.2f)\n\n", data->data[inds].lon[indp1], data->data[inds].lat[indp1]);
*/
  //calculate the side lengths:
  s2=sdist(data->data[inds1].lon[indp], data->data[inds1].lat[indp],
		  data->data[inds1].lon[indp1], data->data[inds1].lat[indp1]);
  s3=sdist(data->data[inds].lon[indp1], data->data[inds].lat[indp1],
		  data->data[inds1].lon[indp1], data->data[inds1].lat[indp1]);

  //use the distances to calculate the weights:
  weight[0]=(d1-d2-s1)*(d1-d4-s4)/s1/s4/4;
  weight[1]=(d2-d1-s1)*(d2-d3-s2)/s1/s2/4;
  weight[2]=(d3-d2-s2)*(d3-d4-s3)/s2/s3/4;
  weight[3]=(d4-d3-s3)*(d4-d1-s4)/s3/s4/4;
  
  //printf("%8.3g %8.3g %8.3g %8.3g\n", weight[0], weight[1], weight[2], weight[3]);

  //can't forget the values:
  val[0]=data->data[inds].bt[indp];
  val[1]=data->data[inds1].bt[indp];
  val[2]=data->data[inds1].bt[indp1];
  val[3]=data->data[inds].bt[indp1];

#ifdef PRINT_FLOC 
  printf("scan indices:  %d %d\n", inds, inds1);
  printf("track indices: %d %d\n", indp, indp1);
#endif

  return;

  printf("distances:\n\n");
  printf("%f\n", sqrt(d1));
  printf("%f\n", sqrt(d2));
  printf("%f\n", sqrt(d3));
  printf("%f\n\n", sqrt(d4));

  printf("side lengths:\n\n");
  printf("%f\n", sqrt(s1));
  printf("%f\n", sqrt(s2));
  printf("%f\n", sqrt(s3));
  printf("%f\n\n", sqrt(s4));

  printf("s11=%f\n", -(d2-d1-s1)/sqrt(s1)/2);
  printf("s12=%f\n", -(d1-d2-s1)/sqrt(s1)/2);

}

amsu_1c_data * amsu_interp_base::get_data_forward(long index) {
  if (amsu_data[index] != NULL) return amsu_data[index];

  //search from zero-point until index for data to delete:
  if (nloaded >= maxload) for (long i=0; i<index; i++) {
    if (amsu_data[i] != NULL) {
      delete_amsu_1c_data(amsu_data[i]);
      amsu_data[i]=NULL;
      nloaded--;
      break;
    }
  }

  //if none found, search backwards from top:
  if (nloaded >= maxload) for (long i=nfiles-1; i>index; i--) {
    if (amsu_data[i] != NULL) {
      delete_amsu_1c_data(amsu_data[i]);
      amsu_data[i]=NULL;
      nloaded--;
      break;
    }
  }


  amsu_data[index]=(*amsu_reader)(filelist[index], SWAP_ENDIAN);
  nloaded++;

  return amsu_data[index];

}

amsu_1c_data *amsu_interp_base::get_data_backward(long index) {
  if (amsu_data[index] != NULL) return amsu_data[index];

  if (nloaded >= maxload) for (long i=nfiles-1; i>index; i--) {
    if (amsu_data[i] != NULL) {
      delete_amsu_1c_data(amsu_data[i]);
      amsu_data[i]=NULL;
      nloaded--;
      break;
    }
  }

  if (nloaded >= maxload) for (long i=0; i<index; i++) {
    if (amsu_data[i] != NULL) {
      delete_amsu_1c_data(amsu_data[i]);
      amsu_data[i]=NULL;
      nloaded--;
      break;
    }
  }

  amsu_data[index]=(*amsu_reader)(filelist[index], SWAP_ENDIAN);
  nloaded++;

  return amsu_data[index];

}

