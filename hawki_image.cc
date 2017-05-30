#include <sstream>
#include "hawki_image.h"


//
//
int HAWKIimage::flatFieldImage() {
  FITS* raw  = this->raw();
  FITS* dark = this->dark();
  FITS* flat = this->flat();
  FITS* bpm  = this->bpm();

  // Create primary HDU based on the primary HDU of the raw image
  //
  const string flatname = getdbimdir() + flatfieldedname() + string(".fits");
  remove(flatname.c_str());  // Delete old file if it already exists
  FITS flatfielded(flatname,(*raw));

  // Create a saturation map
  const string saturname = getdbimdir() + "/saturmap.fits";
  remove(saturname.c_str());
  FITS saturmap(saturname,Write,false);
  
  // Determine exposure time
  //
  /*
    PHDU& rphdu = raw->pHDU();
    float exptime;
    rphdu.readKey("EXPTIME",exptime);
  */

  // Dark subtract and Flatfield one HDU at a time
  // 
  for (int ext=1; ext<=NHDU; ext++) {
    ExtHDU& rawext  = raw->extension(ext);
    ExtHDU& darkext = dark->extension(ext);
    ExtHDU& flatext = flat->extension(ext);
    ExtHDU& bpmext  = bpm->extension(ext);

    // Copy the extension from the raw image to copy the
    // extension header.
    //
    flatfielded.copy(rawext);
    ExtHDU& flatfieldedext = flatfielded.extension(ext);

    valarray<float> fraw, fdark, fflat, fbpm;
    rawext.read(fraw);
    darkext.read(fdark);
    flatext.read(fflat);
    bpmext.read(fbpm);
    
    // All images should have the same size
    //
    long nx(bpmext.axis(0));
    long ny(bpmext.axis(1));

    valarray<unsigned int> satur;
    satur.resize(nx*ny);
    satur = 0.0;

    // Gain in e/s
    //
    float gain = 1.0;
    if (rawext.name() == "CHIP1.INT1")
      gain = 1.705;
    else if (rawext.name() == "CHIP2.INT1")
      gain = 1.870;
    else if (rawext.name() == "CHIP3.INT1")
      gain = 1.735;
    else if (rawext.name() == "CHIP4.INT1") 
      gain = 2.110;

    for (unsigned int n=0; n<fraw.size(); n++) {
      if (fraw[n] >= 60000.0/gain)                // Pixel is saturated
	satur[n] = 1;

      if (fbpm[n] || fflat[n] < 0.05) {           // Check if flat is ~0
	fbpm[n] = 1;                              // Mask pixel as bad
      } else {
	fraw[n] -= fdark[n];                      // Dark subtract
	fraw[n] /= fflat[n];                      // Flat field
      }
    }

    // Gain multiply
    //
    fraw *= gain;

    // Write dark subtracted and flat fielded extension
    //
    vector<long> container(2);
    container[0] = nx;
    container[1] = ny;
    long fpixel(1);
    long nelem = nx * ny;
    flatfieldedext.write(fpixel,nelem,fraw);

    // Write saturation map
    //
    ExtHDU* saturext = saturmap.addImage(rawext.name(),USHORT_IMG,container);
    saturext->write(fpixel,nelem,satur);

    /*
     * We take care of this else where, plus we do not want to
     * mess with the dead.fits file at this stage.
     *
    // There appears to be something going on with the four
    // pixels surrounding each chip (this is probably already
    // masked).
    //
    for (int j=0; j<4; j++) {
      for (int i=0; i<nx; i++) {
	fbpm[nx*j+i] = 1;
	fbpm[nx*(ny-1-j)+i] = 1;
      }
      for (int i=0; i<ny; i++) {
	fbpm[nx*i+j] = 1;
	fbpm[nx*i+nx-1-j] = 1;
      }
    }
    */

    // Write "dead" image.  This is a way of avoiding to overwrite the
    // original bpm image and to make it work with poloka without
    // changing it.
    //
    // deadimageext.write(fpixel,nelem,fbpm);
  }

  // Update header
  //
  PHDU& phdu = flatfielded.pHDU();
  //  phdu.addKey("TOADPIXS",0.106,"Pixel scale arcsec/pixel");
  phdu.addKey("TOADGAIN",1,"Image is in electrons");
  // phdu.addKey("EXPTIME",1,"Image is in e/s");

  flatfielded.destroy();
  
  return 0;
}
