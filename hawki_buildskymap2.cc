/*
 * USAGE
 *
 *   hawki_buildskymap2 <dbimage>
 *
 * ABOUT
 *
 *   Create a skymap from the skysubtracted image by calculating the sky
 *   from a box centered on each pixel.
 *
*/

#include "hawki_image.h"
#include "stats.h"

#include <string>
#include <iostream>
#include <valarray>

using namespace std;


// The minimum number of frames required for trusting a given sky value
//
const int minimused = 3;
const int NHDU = 4;

// Box width = 2*hkw+1
//
const unsigned int hkw = 15;

int main(int argc, char **argv) {
  
  for (int dbi=1; dbi<argc; dbi++) {
    string dir = string(argv[dbi]);
    HAWKIimage dbimage(dir);
    
    cout << dir << endl;

    FITS* f = dbimage.skysubtracted();
    FITS* s = dbimage.skymap();
    FITS* n = dbimage.nimuse();
    FITS* d = dbimage.bpm();
    FITS* o = dbimage.objectmask();

    string skymap2file = dbimage.getdbimdir() +
      string("/skymap2.fits");
    remove(skymap2file.c_str());  // Delete old file if it already exists
    FITS skymap2fits(skymap2file,Write,false);
    
    for (int ext=1; ext<=NHDU; ext++) {
      ExtHDU& fHDU = f->extension(ext);
      ExtHDU& sHDU = s->extension(ext);
      ExtHDU& nHDU = n->extension(ext);
      ExtHDU& dHDU = d->extension(ext);
      ExtHDU& oHDU = o->extension(ext);

      string extname;
      extname.assign(fHDU.name());

      vector<long> extAx(2);
      extAx[0] = fHDU.axis(0);
      extAx[1] = fHDU.axis(1);

      long fpixel(1);
      long npix = extAx[0] * extAx[1];

      valarray<float> image, skymap, nimuse, omask, dead;
      fHDU.read(image);
      sHDU.read(skymap);
      nHDU.read(nimuse);
      dHDU.read(dead);
      oHDU.read(omask);
      
      //  Update mask
      //
      for (unsigned int i=0; i<image.size(); i++) {
	// Mask all pixels where the sky value is not based on at least
	// minimused frames, and the pixel value is less than min value
	//
	if ( nimuse[i] < minimused || dead[i] == 1 ) omask[i] = 1;
      }


      valarray<float> skymap2(npix);

      // Calcualte skymap
      //
      for (unsigned int ii=0; ii<extAx[0]; ii++) {

	for (unsigned int jj=0; jj<extAx[1]; jj++) {

	  unsigned int i = ii*extAx[0] + jj;
	  
	  vector<float> sky;
	  if (dead[i] == 0) {
	    for (unsigned int ki=ii-hkw; ki<=ii+hkw; ki++) {
	      for (unsigned int kj=jj-hkw; kj<=jj+hkw; kj++) {

		unsigned int j = ki*extAx[0] + kj;
		if (ki > 0 && kj > 0 && ki < extAx[0] && kj < extAx[1] &&
		    omask[j] == 0. )
		  sky.push_back((float) image[j]);
	      }
	    }
	  }
	    
	  
	  // Update smoothed sky map with the median value
	  //
	  if (sky.size() > 3 && dead[i] == 0) {
	    sort(&sky[0],&sky[0] + sky.size());
	    if (sky.size() % 2)
	      skymap2[i] = sky[(sky.size()-1)/2];
	    else
	      skymap2[i] = 0.5*(sky[sky.size()/2-1] + sky[sky.size()/2]);

	  } else {
	    skymap2[i] = 0.;
	  }
	  
	}
      }

      // Write smoothed sky frame
      //
      ExtHDU* skymap2ext = skymap2fits.addImage(extname,FLOAT_IMG,extAx);
      skymap2ext->write(fpixel,npix,skymap2);
    }

    skymap2fits.destroy();
  }
}

// EOF
//
