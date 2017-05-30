/*
 * USAGE
 *
 *   hawki_skysubtract <dbimage>
 *
 * ABOUT
 *
 *   Subtract the skymap from hawki_buildskymap from the flatfielded image
 *   (from hawki_flatfield) and produce a sky subtracted frame and a mask
 *   which is the combination of the bad pixel map and the areas where a
 *   good sky estimate could not be obtained.
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

int main(int argc, char **argv) {

  cerr << "WARNING: WILL MASK ALL VALUES BELOW " << minvalue << 
    " BEFORE SKY SUBTRACTING!" << endl;
  
  for (int dbi=1; dbi<argc; dbi++) {
    string dir = string(argv[dbi]);
    HAWKIimage dbimage(dir);
    
    cout << dir << endl;

    FITS* f = dbimage.flatfielded();
    FITS* s = dbimage.skymap();
    FITS* n = dbimage.nimuse();
    FITS* d = dbimage.bpm();
    FITS* o = dbimage.objectmask();

    string skysubfile = dbimage.getdbimdir() +
      string("/skysubtracted.fits");
    remove(skysubfile.c_str());  // Delete old file if it already exists
    FITS skysubfits(skysubfile,*f);

    string badfile = dbimage.getdbimdir() +
      string("/mask.fits");
    remove(badfile.c_str());
    FITS badfits(badfile,*d);
    
    string weightfile = dbimage.getdbimdir() + 
      string("/skysubtracted.weight.fits");
    remove(weightfile.c_str());
    FITS weightfits(weightfile,*f);

    for (int ext=1; ext<=NHDU; ext++) {
      ExtHDU& fHDU = f->extension(ext);
      ExtHDU& sHDU = s->extension(ext);
      ExtHDU& nHDU = n->extension(ext);
      ExtHDU& dHDU = d->extension(ext);

      skysubfits.copy(fHDU);
      weightfits.copy(fHDU);
      badfits.copy(dHDU);

      valarray<float> image, skymap, nimuse, omask, dead;
      fHDU.read(image);
      sHDU.read(skymap);
      nHDU.read(nimuse);
      dHDU.read(dead);      
      if (o != NULL) {
	ExtHDU& oHDU = o->extension(ext);
	oHDU.read(omask);
      } else {
	omask.resize(dead.size());
	omask = 0;
      }

      // Skysubtract and update bad pixel map
      //
      for (unsigned int i=0; i<image.size(); i++) {
	// Mask all pixels where the sky value is not based on at least
	// minimused frames, and the pixel value is less than min value
	//
	if (nimuse[i] < minimused || image[i] < minvalue) dead[i] = 1;

	if (dead[i] == 0)
	  image[i] -= skymap[i];
	else
	  image[i] = 0.0;
      }

      // Set up a weight map to the inverted dead map for now, we
      // will update the non-masked pixels with the inverse variance
      // below.
      //
      valarray<float> weight(dead.size());
      weight += 1.0;
      weight -= dead;

      // Write output files
      //
      vector<long> extAx(2);
      extAx[0] = fHDU.axis(0);
      extAx[1] = fHDU.axis(1);
      long fpixel(1);
      long npix = extAx[0] * extAx[1];

      // Write the updated dead pixel map
      //
      ExtHDU& badext = badfits.extension(ext);
      badext.write(fpixel,npix,dead);

      // Get the original sky value
      //
      sHDU.readAllKeys();
      Keyword& osky = sHDU.keyWord("ORIGSKY");
      float origsky;
      osky.value(origsky);

      // Create some sky statistics, to make sure that the
      // sky subtraction was successful.
      //
      dead += omask;   // Updated with object mask, we only want sky pixels!
      float median, mean, std_dev;
      double nelem = image.size();
      clippedStats<3,5>(image,dead,median,mean,std_dev);
      printf("\tEXT: %d, origsky = %.2f\n", ext, origsky);
      printf("\t\t(mean,rms,median,nelem,rms/sqrt(nelem)) = (%.3f,%.3f,%.3f,%d,%.3f)\n",
	     mean,std_dev,median,(int) nelem,std_dev/sqrt(nelem));
      

      // Write sky subtracted frame
      //
      ExtHDU& skysubext = skysubfits.extension(ext);

      skysubext.addKey("SKYLEV",origsky,"Original sky level");
      skysubext.addKey("SKYSIGEX",std_dev,"Original sky sigma");
      skysubext.addKey("SKYSIGTH",sqrt(origsky),"Theoretical sigma");

      skysubext.addKey("BACKLEV",0,"");
      skysubext.addKey("BACK_SUB",true,"");
      skysubext.write(fpixel,npix,image);

      // Write weight frame after we set all non-dead pixels to the
      // inverted variance.
      //
      weight *= 1.0/(std_dev*std_dev);
      ExtHDU& weightext = weightfits.extension(ext);
      weightext.write(fpixel,npix,weight);
    }

    // Add the product type for the ESO reduction pipeline (should
    // use setValue instead).
    //
    PHDU& phdu = skysubfits.pHDU();
    phdu.addKey("HIERARCH ESO PRO CATG","BKG_SUBTRACTED","");
    //    Keyword& procatg = phdu.keyWord("HIERARCH ESO PRO CATG");
    //    procatg.setValue(string("BKG_SUBTRACTED"));

    badfits.destroy();
    skysubfits.destroy();
    weightfits.destroy();
  }
}

// EOF
//
