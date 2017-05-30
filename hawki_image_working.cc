#include <sstream>
#include "hawki_image.h"

int HAWKIimage::flatfieldimage() {
  auto_ptr<FITS> raw = readRaw();
  auto_ptr<FITS> dark = readDark();
  auto_ptr<FITS> flat = readFlat();

  PHDU& praw = raw->pHDU();

  // Create empty primary HDU
  //
  const string filename = getdbimdir() + "/flatfielded.fits";
  auto_ptr<FITS> flatfielded(0);
  try {
    long naxes[2] = { 0, 0 };
    remove(filename);              // Delete old file if it already exists
    flatfielded.reset(new FITS(filename,USHORT_IMG,0,naxes));
  } catch (FITS::CantCreate) {
    cerr << "ERROR: failed to create :" << filename << endl;
  }
  
  // Dark subtract and Flatfield one HDU at a time
  // 
  for (int ext=1; ext<=NHDU; ext++) {
    ExtHDU& rawext = raw->extension(ext);
    ExtHDU& darkext = dark->extension(ext);
    ExtHDU& flatext = flat->extension(ext);

    valarray<float> fraw, fdark, fflat;
    rawext.read(fraw);
    darkext.read(fdark);
    flatext.read(fflat);
    
    fraw -= fdark;
    fraw /= fflat;

    rawext.readAllKeys();

    stringstream extname;               // Convert to string
    extname << ext;

    vector<long> container(2,rawext.axis(1));
    ExtHDU* flatfieldedext =
      flatfielded->addImage(extname.str(),FLOAT_IMG,container);
    long fpixel(1);
    long nelem = rawext.axis(0) * rawext.axis(1);
    //    flatfieldedext->copyAllKeys(&praw);
    flatfieldedext->write(fpixel,nelem,fraw);
  }

  //  PHDU& phdu = flatfielded->pHDU();
  //  (&phdu)->copyAllKeys(&(raw->pHDU()));
  //  phdu.write(long(1));

  return 0;
}
