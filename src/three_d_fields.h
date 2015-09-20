#define RDRYAIR 287
#define RWATERVAPOUR 461
#define GACC 9.8

float *** allocate_3d_field(long nlev, long nlat, long nlon);
void delete_3d_field(float *** field, long nlev, long nlat, long nlon);
float *** read_grads_bin(char *filename, long nlev, long nlat, long nlon);
int read_grads_bin(char *filename, long nlev, long nlat, long nlon, 
		float ***data);
float *** calc_height(float *** pp, float *** tt, float *** vmr, 
		long nlev, long nlat, long nlon);


