#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <getopt.h>

#include "peteys_tmpl_lib.h"
#include "linked.cc"
#include "read_gome.h"

#define MAX_LL 1000

#define NORMAL 0
#define RAW 1

#define HEADER_LENGTH 200

//#define ASCII_OUT

//describes scan line with at least one measurement in the bin:
struct amsu_bin_record {
  time_class date;
  long file_index;
  long scan_index;
};


int read_global_bins (char *fname, long *&nlon, long &nlat, float &offset, float &fudge_factor) {
  FILE *fs;
  float lat;
  float checklat;
  long err=0;

  fs=fopen(fname, "r");
  fscanf(fs, "%ld", &nlat);
  nlon=new long [nlat+2];
  fscanf(fs, "%f", &fudge_factor);
  fscanf(fs, "%f %ld", &lat, nlon+1);
  offset=90+lat;
  for (long i=1; i<nlat; i++) {
    checklat=i*(180-offset*2)/nlat-90+offset;
    fscanf(fs, "%f %d", &lat, nlon+i+1);
    if (lat!=checklat) {
      fprintf(stderr, "Scanned and calculated latitudes do not agree: %f vs. %f\n", lat, checklat);
      err=1;
    }
  }

  nlon[0]=1;
  nlon[nlat+1]=1;

  fclose(fs);

  return err;

}

void set_default_bins2(long *&nlon, long &nlat, float &offset, float &fudge_factor) {
  nlat=17;
  offset=5;
  fudge_factor=0.05;

  nlon=new long[nlat+2];
  nlon[0]=1;
  nlon[1]=9;
  nlon[2]=15;
  nlon[3]=21;
  nlon[4]=26;
  nlon[5]=30;
  nlon[6]=33;
  nlon[7]=35;
  nlon[8]=36;
  nlon[9]=36;
  nlon[10]=36;
  nlon[11]=35;
  nlon[12]=33;
  nlon[13]=30;
  nlon[14]=26;
  nlon[15]=21;
  nlon[16]=15;
  nlon[17]=9;
  nlon[18]=1;

}

void set_default_bins(long *&nlon, long &nlat) {
  float mindlon=50;
  float maxdlon=111.111111;
  float lat;
  float dlon;
  
  nlat=180;
  nlon=new long[nlat];
  nlon[0]=1;
  nlon[nlat-1]=1;
  for (long i=1; i<nlat/2; i++) {
    dlon=maxdlon+i*(mindlon-maxdlon)/89;
    lat=M_PI*i/nlat;
    nlon[nlat/2-i]=(long) (40000.*cos(lat)/dlon);
    if (nlon[nlat/2-i] > 360) nlon[nlat/2-i]=360;
    nlon[nlat/2+i-1]=nlon[nlat/2-i];
  }

  for (long i=0; i<nlat; i++) {
    lat=i*180/nlat-90;
    printf("%8.2f %5d\n", lat, nlon[i]);
  }
}

int main(int argc, char ** argv) {
  char *initfile;
  char *indexfile;
  char *binfile;
  char c;

  char checkfile[MAX_LL];

  int swap_flag;
  int amsu_read_type;

  int4_t (*amsu_head_reader) (char *, time_class &, time_class &, int4_t &, int);
  gome_data * (*amsu_reader) (char *, int);

  char *basepath;
  long baselen;
  char **filelist;
  char filename[MAX_LL];
  long nfiles;
  time_class *tstart;
  time_class *tend;
  long *nlon, nlat;

  float offset;
  float fudge_factor;

  linked_list<amsu_bin_record *> *** bins;
  amsu_bin_record * bin_rec;

  FILE *fs, *fs1;
  char line[MAX_LL];
  long fp_amsustart;		//start of list of amsu files
  long ll;			//line length
  char ** flist1;		//for sorting file-list
  time_class *t1, *t2;
  long *sind;			//sorting indices

  //throw-away:
  long lonind, latind;
  long *binind;
  gome_data *data;
  int4_t nchan;
  char ts1[30], ts2[30];

  long nent;		//number of entries in the bin
  amsu_bin_record **bin_contents;
  amsu_bin_record **bin_contents2;
  time_class tlast;	//for removing duplicate entries

  int *dup;
  long ndup;
  long nrec;
  long ndup_total;

  float lon, lat;

  char header[HEADER_LENGTH];

  float lon_offset;

  //set defaults before parsing command line:
  binfile=NULL;
  basepath=NULL;

  swap_flag=0;

  amsu_read_type=RAW;

  lon_offset=180;
  
  //parse the command line arguments:
  while ((c = getopt (argc, argv, "0Rsf:p:")) != -1) {
    switch (c) {
        case ('R'):
	    amsu_read_type=RAW;
            break;
        case ('s'):
            swap_flag=1;
            break;
	case ('f'):
	    binfile=new char[strlen(optarg)+1];
	    strcpy(binfile, optarg);
	    break;
	case ('p'):
	    basepath=new char[strlen(optarg)+1];
	    strcpy(basepath, optarg);
	    break;
	case ('0'):
            lon_offset=0;
	    break;
        case ('?'):
            fprintf(stderr, "Unknown option: %c --ignored\n", optopt);
            break;
        default:
            fprintf(stderr, "Error parsing command line\n");
            exit(2);
    }
  }

  if (argc-optind < 2) {
    printf("usage: index_gome [-R] [-s] [-f binfile] [-p basepath]\n");
    printf("                  initfile indexfile\n");
    printf("\n");
    printf("arguments:\n");
    printf("  initfile     = list of AMSU data files to be indexed\n");
    printf("  indexfile    = output index file\n");
    printf("\n");
    printf("options\n");
    printf("  -R           = use \"raw\" read routines\n");
    printf("  -s           = byte-swap the files\n");
    printf("  -0           = 0 to 360 degree longitude convention\n");
    printf("                 (default is -180 to 180 degrees)\n");
    printf("  -f binfile   = set bins from file\n");
    printf("  -p basepath  = base directory of AMSU l1c files\n");
    printf("\n");
    exit(1);
  }

  if (amsu_read_type == RAW) {
    amsu_head_reader=&rawread_gome_head;
    amsu_reader=&rawread_gome;
  } else {
//    amsu_head_reader=&read_amsu_1c_head;
//    amsu_reader=&read_amsu_1c;
  }

  if (binfile != NULL) {
    read_global_bins(binfile, nlon, nlat, offset, fudge_factor);
  } else {
    set_default_bins2(nlon, nlat, offset, fudge_factor);
  }

  if (basepath == NULL) {
    baselen=0;
    strcpy(filename, "");
  } else {
    baselen=strlen(basepath);
    strcpy(filename, basepath);
  }

  initfile=argv[optind];
  indexfile=argv[optind+1];

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
  for (long i=0; i<nfiles; i++) {
    //read the name of the file:
    fgets(line, MAX_LL, fs);
    ll=strlen(line);
    line[ll-1]=0;
    filelist[i]=new char[ll];
    strcpy(filelist[i], line);

    //read start and end dates:
    strcpy(filename+baselen, filelist[i]);
    (*amsu_head_reader)(filename, tstart[i], tend[i], nchan, swap_flag);
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

  //create the bins:
  bins=new linked_list<amsu_bin_record *> ** [nlat+2];
  for (long i=0; i<nlat+2; i++) {
    bins[i]=new linked_list<amsu_bin_record *> * [nlon[i]];
    for (long j=0; j<nlon[i]; j++) bins[i][j]=new linked_list<amsu_bin_record *>;
  }


  //bin each of the data files one by one:
  for (long ifile=0; ifile<nfiles; ifile++) {
    strcpy(filename+baselen, filelist[ifile]);
    data=(*amsu_reader) (filename, swap_flag);
    binind=new long [9];
    if (data == NULL) {
      fprintf(stderr, "Error: cannot read input file %s\n", filelist[ifile]);
      exit(5);
    }
    for (long i=0; i<data->npix; i++) {
      //we want to bin each scan line as few times as possible:
      latind=(long) ((data->data[i].lat+90.-offset)/(180.-2*offset)*nlat+1);
      lonind=(long) ((data->data[i].lon+lon_offset)/360.*nlon[latind]);
      //printf("%f %f %d %d %f %d\n", data->data[i].lon, data->data[i].lat, lonind, latind, offset, nlat);
      if (lonind == nlon[latind]) lonind=nlon[latind]-1;
      binind[0]=lonind*(nlat+2)+latind;
	
      lonind=(long) ((data->data[i].lon+lon_offset)/360.*nlon[latind]-fudge_factor);
      if (lonind < 0) lonind=0;
      if (lonind == nlon[latind]) lonind=nlon[latind]-1;
      binind[7]=lonind*(nlat+2)+latind;

      lonind=(long) ((data->data[i].lon+lon_offset)/360.*nlon[latind]+fudge_factor);
      if (lonind >= nlon[latind]) lonind=nlon[latind]-1;
      binind[8]=lonind*(nlat+2)+latind;

      latind=(long) ((data->data[i].lat+90.-offset)/(180.-2*offset)*nlat-fudge_factor+1);
      lonind=(long) ((data->data[i].lon+lon_offset)/360.*nlon[latind]);
      if (lonind == nlon[latind]) lonind=nlon[latind]-1;
      binind[1]=lonind*(nlat+2)+latind;
	
      lonind=(long) ((data->data[i].lon+lon_offset)/360.*nlon[latind]-fudge_factor);
      if (lonind < 0) lonind=0;
      if (lonind == nlon[latind]) lonind=nlon[latind]-1;
      binind[2]=lonind*(nlat+2)+latind;

      lonind=(long) ((data->data[i].lon+lon_offset)/360.*nlon[latind]+fudge_factor);
      if (lonind >= nlon[latind]) lonind=nlon[latind]-1;
      binind[3]=lonind*(nlat+2)+latind;

      latind=(long) ((data->data[i].lat+90.-offset)/(180.-2*offset)*nlat+fudge_factor+1);
      lonind=(long) ((data->data[i].lon+lon_offset)/360.*nlon[latind]);
      if (lonind == nlon[latind]) lonind=nlon[latind]-1;
      binind[4]=lonind*(nlat+2)+latind;

      lonind=(long) ((data->data[i].lon+lon_offset)/360.*nlon[latind]-fudge_factor);
      if (lonind < 0) lonind=0;
      if (lonind == nlon[latind]) lonind=nlon[latind]-1;
      binind[5]=lonind*(nlat+2)+latind;

      lonind=(long) ((data->data[i].lon+lon_offset)/360.*nlon[latind]+fudge_factor);
      if (lonind >= nlon[latind]) lonind=nlon[latind]-1;
      binind[6]=lonind*(nlat+2)+latind;

      //sort the bin indices so we can throw away duplicates:
      heapsort_inplace(binind, 9);

      //add the new index entry:
      bin_rec=new amsu_bin_record;
      bin_rec->date=data->data[i].date;
      bin_rec->file_index=ifile;
      bin_rec->scan_index=i;
      latind=binind[0] % (nlat+2);
      lonind=binind[0]/(nlat+2);
      bins[latind][lonind]->add(bin_rec);

      for (long j=1; j<9; j++) {
        if (binind[j] != binind[j-1]) {
          //wastes space, should really do this with objects...
          bin_rec=new amsu_bin_record;
          bin_rec->date=data->data[i].date;
          bin_rec->file_index=ifile;
          bin_rec->scan_index=i;
          latind=binind[j] % (nlat+2);
          lonind=binind[j]/(nlat+2);
          bins[latind][lonind]->add(bin_rec);
	}
      }

    }
    delete_gome_data(data);
    delete [] binind;
  }

  //write out the index to a file:
  //fs=stdout;

  fs1=fopen(indexfile, "w");
  sprintf(header, "Index of data points in AMSU files in variable size bins");
  fwrite(header, sizeof(char), HEADER_LENGTH-sizeof(float), fs1);
  //try to keep it backwards compatible:
  lon_offset=lon_offset-180;
  fwrite(&lon_offset, sizeof(float), 1, fs1);
  lon_offset=lon_offset+180;

#ifdef ASCII_OUT
  fs=fopen("check_index.txt", "w");
  fprintf(fs, "%s\n", header);
  if (amsu_read_type == RAW) {
    fprintf(fs, "%4d RAW files\n", nfiles);
  } else {
    fprintf(fs, "%4d NORMAL files\n", nfiles);
  }
  fprintf(fs, "%s\n", basepath);
#endif

  fwrite(&amsu_read_type, sizeof(amsu_read_type), 1, fs1);
  fwrite(&swap_flag, sizeof(swap_flag), 1, fs1);
  fwrite(basepath, sizeof(char), baselen+1, fs1);
  fwrite(&nfiles, sizeof(nfiles), 1, fs1);
  for (long i=0; i<nfiles; i++) {

#ifdef ASCII_OUT	  
    fprintf(fs, "%s\n", filelist[i]);
#endif

    fwrite(filelist[i], sizeof(char), strlen(filelist[i])+1, fs1);
  }

  fwrite(&nlat, sizeof(long), 1, fs1);
  fwrite(nlon, sizeof(long), nlat+2, fs1);
  ndup_total=0;
  for (long i=0; i<nlat+2; i++) {
    lat=(i-1)*(180.-2*offset)/nlat-90.+offset;
    if (i == 0) lat=-90;
    printf("%d\n", nlon[i]);
    for (long j=0; j<nlon[i]; j++) {
      lon=j*360./nlon[i]-lon_offset;

#ifdef ASCII_OUT
      sprintf(checkfile, "check_bins/bin(%7.2f,%7.2f).txt", lon, lat);
      fs=fopen(checkfile, "w");
#endif

      bin_contents2=bins[i][j]->make_array(nent);
      printf("(%d, %d); (%f, %f): %d\n", i, j, lon, lat, nent);
      if (nent == 0) {
        fprintf(stderr, "Warning: empty bin found: (%d, %d); (%f, %f)\n", i, j, lon, lat);
#ifdef ASCII_OUT	  
        //fprintf(fs, "(%6d, %6d); (%8.2f, %8.2f): %d\n", i, j, lon, lat, nent);
#endif
	fwrite(&i, sizeof(long), 1, fs1);
	fwrite(&j, sizeof(long), 1, fs1);
	fwrite(&lon, sizeof(float), 1, fs1);
	fwrite(&lat, sizeof(float), 1, fs1);
	fwrite(&nent, sizeof(nent), 1, fs1);

        delete bins[i][j];
#ifdef ASCII_OUT	 
        fclose(fs);
#endif 
	continue;
      }

      //sort by date:
      t1=new time_class[nent];
      for (long k=0; k<nent; k++) t1[k]=bin_contents2[k]->date;
      sind=heapsort(t1, nent);
      bin_contents=new amsu_bin_record * [nent];
      for (long k=0; k<nent; k++) bin_contents[k]=bin_contents2[sind[k]];

      //search for duplicates:
      dup=new int[nent];
      dup[0]=0;
      ndup=0;
      tlast=bin_contents[0]->date;
      for (long k=1; k<nent; k++) {
        dup[k]=0;
        if (bin_contents[k]->date==tlast) {
	  fprintf(stderr, "Duplicate record found and removed\n");
          printf("(%4d, %4d): %d\n", i-180, j-90, nent);
          printf("%29s %4d %10d\n", ts1, bin_contents[k-1]->file_index, bin_contents[k-1]->scan_index);
          printf("%29s %4d %10d\n", ts1, bin_contents[k]->file_index, bin_contents[k]->scan_index);
	  dup[k]=1;
	  ndup++;
	}
        tlast=bin_contents[k]->date;
      }

#ifdef ASCII_OUT	  
      //fprintf(fs, "(%6d, %6d); (%8.2f, %8.2f): %d\n", i, j, lon, lat, nent-ndup);
#endif
      nrec=nent-ndup;
      ndup_total+=ndup;
      fwrite(&i, sizeof(long), 1, fs1);
      fwrite(&j, sizeof(long), 1, fs1);
      fwrite(&lon, sizeof(float), 1, fs1);
      fwrite(&lat, sizeof(float), 1, fs1);
      fwrite(&nrec, sizeof(long), 1, fs1);
      for (long k=0; k<nent; k++) {
        if (dup[k]==1) continue;
#ifdef ASCII_OUT	  
        bin_contents[k]->date.write_string(ts1);
        fprintf(fs, "%29s %4d %10d\n", ts1, bin_contents[k]->file_index, bin_contents[k]->scan_index);
#endif
	fwrite(&bin_contents[k]->date, sizeof(time_class), 1, fs1);
      }
      for (long k=0; k<nent; k++) {
        if (dup[k]==1) continue;
	fwrite(&bin_contents[k]->file_index, sizeof(long), 1, fs1);
      }
      for (long k=0; k<nent; k++) {
        if (dup[k]==1) continue;
	fwrite(&bin_contents[k]->scan_index, sizeof(long), 1, fs1);
      }
      //fflush(fs1);

      delete [] t1;
      for (long k=0; k<nent; k++) delete bin_contents[k];
      delete [] bin_contents;
      delete [] bin_contents2;
      delete [] sind;
      delete [] dup;

      delete bins[i][j];
#ifdef ASCII_OUT	 
      fclose(fs);
#endif 
    }
    delete [] bins[i];
  }
  fprintf(stderr, "%d duplicates found\n", ndup_total);  
  fclose(fs1);
#ifdef ASCII_OUT	 
  //fclose(fs);
#endif 

  delete [] bins;  
  for (long i=0; i<nfiles; i++) delete [] filelist[i];
  delete [] filelist;
  delete [] tstart;
  delete [] tend;
  delete [] nlon;

  if (binfile != NULL) delete [] binfile;
  if (basepath != NULL) delete [] basepath;

}


