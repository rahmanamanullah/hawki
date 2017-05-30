#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

/* Timing function */
double gettime(void) {
  struct timeval tv;
  gettimeofday(&tv,NULL);
  return tv.tv_sec + 1e-6*tv.tv_usec;
}

/* Running mean */
void runningmean(int nim, int nrej, float *imagpx, float *maskpx,
		 float *sky, float *rms, int *nimused) {
  int i,ii,jj,n,n1,n2;
  int mval, skysize;
  float foo;
  float pval, *skyvec;
  int inserted;

  skyvec = (float *) malloc(sizeof(float) * nim);

  skysize = 0;
  for (i = 0; i < nim; i++) {  // Loop over image stack
    pval = imagpx[i];          // Data value
    mval = maskpx[i];          // Mask value 

    //    printf("%.2f %d\n", pval, mval);
	  
    // Insert into a sorted list unless the pixel is masked
    //
    if (mval == 0) {                                  // Not masked
      inserted = 0;
      ii = 0;
      while (inserted == 0 && skysize > 0 && ii < skysize) {
	if (pval < skyvec[ii]) {
	  for (jj=skysize-1; jj>=ii; jj--)            // Shift rest of vector
	    skyvec[jj+1] = skyvec[jj];
	  skyvec[ii] = pval;                          // Insert value
	  inserted = 1;
	}
	ii++;
      }
      if (inserted == 0) skyvec[skysize] = pval;          // Insert at the end
      skysize++;
    }
  }

  
  *nimused = skysize - 2*nrej;   // Number of images used for this pixel
  if (*nimused < 0) *nimused = 0;
  
  // Are there enough values in the list to continue?
  //
  if (*nimused > 1) {
    n1=nrej;
    n2=skysize-nrej;
    
    // Calculate mean and RMS
    //
    *sky = 0.0;
    *rms = 0.0;
    for (n=n1; n<n2; n++) {
      foo = skyvec[n];
      *sky += foo;
      *rms += foo * foo;
    }
    *sky /= n2-n1;
    *rms /= n2-n1;
    *rms = sqrt(*rms-(*sky)*(*sky));
  } else {
    *sky = 0.0;
    *rms = 0.0;
  }

  //  printf("%.3f %.3f %d\n", *sky, *rms, *nimused);

  free(skyvec);
}


// Run runningmean for all pixels
//
void nrunningmean(int npix, int nim, int nrej, float **cimagpx, float **maskpx, float *sky, float *rms, int *nimused) {
  int j, nimval;
  float skyval, rmsval;

  for (j=0; j<npix; j++) {         // Loop over pixels
    runningmean(nim, nrej, cimagpx[j], maskpx[j], &skyval, &rmsval, &nimval);
    sky[j] = skyval;
    rms[j] = rmsval;
    nimused[j] = nimval;
  }
}
