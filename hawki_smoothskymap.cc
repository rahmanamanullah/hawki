/*
 * USAGE
 *
 *   hawki_smoothskymap <dbimage>
 *
 * ABOUT
 *
 *
*/

#include "hawki_image.h"
#include "stats.h"

#include <string>
#include <iostream>
#include <valarray>

using namespace std;


// All values below this in the flat fielded frame will be masked!  Typically,
// the sky is at several thousand counts.
//
const float minvalue = 50.0;


// The minimum number of frames required for trusting a given sky value
//
const int minimused = 3;
const int NHDU = 4;

// Kernel width = 2*hkw+1
//
const unsigned int hkw = 3;

int main(int argc, char **argv) {

  cerr << "WARNING: WILL MASK ALL VALUES BELOW " << minvalue << "!" << endl;
  
  for (int dbi=1; dbi<argc; dbi++) {
    string dir = string(argv[dbi]);
    HAWKIimage dbimage(dir);
    
    cout << dir << endl;

    FITS* f = dbimage.flatfielded();
    FITS* s = dbimage.skymap();
    FITS* n = dbimage.nimuse();
    FITS* d = dbimage.bpm();

    string smoothskyfile = dbimage.getdbimdir() +
      string("/skymap_smooth.fits");
    remove(smoothskyfile.c_str());  // Delete old file if it already exists
    FITS smoothskyfits(smoothskyfile,Write,false);
    
    for (int ext=1; ext<=NHDU; ext++) {
      ExtHDU& fHDU = f->extension(ext);
      ExtHDU& sHDU = s->extension(ext);
      ExtHDU& nHDU = n->extension(ext);
      ExtHDU& dHDU = d->extension(ext);

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
      
      //  Update mask
      //
      for (unsigned int i=0; i<image.size(); i++) {
	// Mask all pixels where the sky value is not based on at least
	// minimused frames, and the pixel value is less than min value
	//
	//	if (nimuse[i] < minimused || image[i] < minvalue) dead[i] = 1;
      }

      valarray<float> smoothsky(npix);

      // Smoothen skymap and update bad pixel map
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
		    dead[j] == 0.)
		  sky.push_back((float) skymap[j]);
	      }
	    }
	  }

	  // Update smoothed sky map with the median value
	  //
	  if (sky.size() > 3 && dead[i] == 0.) {
	    
	    // Calculate median
	    //
	    sort(&sky[0],&sky[0] + sky.size());
	    if (sky.size() % 2)
	      smoothsky[i] = sky[(sky.size()-1)/2];
	    else
	      smoothsky[i] = 0.5*(sky[sky.size()/2-1] + sky[sky.size()/2]);

	  } else {
	    smoothsky[i] = skymap[i];
	  }
	  
	}
      }

      // Write smoothed sky frame
      //
      ExtHDU* smoothskyext = smoothskyfits.addImage(extname,FLOAT_IMG,extAx);
      /*
      Keyword& osky = sHDU.keyWord("ORIGSKY");
      float origsky;
      osky.value(origsky);
      smoothskyext->addKey("ORIGSKY",origsky,"Median sky");
      Keyword& rsky = sHDU.keyWord("SKYRATI");
      float skyratio;      
      rsky.value(skyratio);
      smoothskyext->addKey("SKYRATI",skyratio,"Ratio between image and sky map");
            */
      smoothskyext->write(fpixel,npix,smoothsky);
    }

    smoothskyfits.destroy();
  }
}

// EOF
//
