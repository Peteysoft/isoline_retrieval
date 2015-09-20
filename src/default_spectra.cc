#include <math.h>

#define LAMBDA0 290.
#define LAMBDAN 790.
#define NLAMBDA 1000

//generates a single spectra in a geometric series:
float *default_spectra(float lambda0, float lambdan, long n) {
  float *spec;

  spec=new float[n+1];

  for (long i=0; i<n; i++) {
    spec[i]=lambda0*pow(lambdan/lambda0, (float) i/(float) n);
  }

  return spec;

}
  
//generates a single spectra in a geometric series:
float *default_spectra_centre(float lambda0, float lambdan, long n) {
  float *spec;

  spec=new float[n];

  for (long i=0; i<n; i++) {
    spec[i]=lambda0*pow(lambdan/lambda0, (i+0.5)/(float) n);
  }

  return spec;

}
  
