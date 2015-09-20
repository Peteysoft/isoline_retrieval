#include "read_amsu.h"

void write_amsu_1c(amsu_1c_data *data, 
		time_class t1, 
		time_class t2, 
		int type)
{


  long *rec;
  long reclen;
  long toffset, lloffset, btoffset;

  if (type==AMSU_A) {
    reclen=768;
    toffset=3;
    lloffset=25;
    btoffset=208;
  } else if (type==AMSU_B) {
    reclen=1152;
    toffset=3;
    lloffset=14;
    btoffset=557;
  } else if (type==AMSU_CLASS_RET) {
    reclen=
