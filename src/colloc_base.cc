#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "peteys_tmpl_lib.h"
#include "colloc_base.h"

#include "global_metric.h"

#define MAX_LL 200
#define SWAP_ENDIAN 1

#define COMMENT_LENGTH 200

//#define CHECK_INDEX_FILE 

//#define PRINT_FLOC

gome_data * colloc_base::get_data_forward(long index) {

  if (amsu_data[index] != NULL) return amsu_data[index];

  //search from zero-point until index for data to delete:
  if (nloaded >= maxload) for (long i=0; i<index; i++) {
    if (amsu_data[i] != NULL) {
      delete_gome_data(amsu_data[i]);
      amsu_data[i]=NULL;
      nloaded--;
      break;
    }
  }

  //if none found, search backwards from top:
  if (nloaded >= maxload) for (long i=nfiles-1; i>index; i--) {
    if (amsu_data[i] != NULL) {
      delete_gome_data(amsu_data[i]);
      amsu_data[i]=NULL;
      nloaded--;
      break;
    }
  }

  //tack on the search path to the file name:
  strcpy(fname, basepath);
  strcat(fname, filelist[index]);
  amsu_data[index]=(*amsu_reader)(fname, SWAP_ENDIAN);
  nloaded++;

  return amsu_data[index];

}

gome_data *colloc_base::get_data_backward(long index) {
  if (amsu_data[index] != NULL) return amsu_data[index];

  if (nloaded >= maxload) for (long i=nfiles-1; i>index; i--) {
    if (amsu_data[i] != NULL) {
      delete_gome_data(amsu_data[i]);
      amsu_data[i]=NULL;
      nloaded--;
      break;
    }
  }

  if (nloaded >= maxload) for (long i=0; i<index; i++) {
    if (amsu_data[i] != NULL) {
      delete_gome_data(amsu_data[i]);
      amsu_data[i]=NULL;
      nloaded--;
      break;
    }
  }

  strcpy(fname, basepath);
  strcat(fname, filelist[index]);
  amsu_data[index]=(*amsu_reader)(fname, SWAP_ENDIAN);
  nloaded++;

  return amsu_data[index];

}

