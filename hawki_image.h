#ifndef HAWKI_IMAGE_H
#define HAWKI_IMAGE_H

#include "dbimage.h"

using namespace std;

class HAWKIimage : public DbImage {
 private:
  const static int NHDU = 4;
  FITS *objectmaskpt, *skymappt, *skysubtpt, *nimusept;

 public :
  int idx;
  HAWKIimage() {
    HAWKIimage(string("0"));
    objectmaskpt = NULL;
    skymappt = NULL;
    nimusept = NULL;
    skysubtpt = NULL;
  };
  HAWKIimage(string dir) : DbImage(dir), idx(0) {
    objectmaskpt = NULL;
    skymappt = NULL;
    nimusept = NULL;
    skysubtpt = NULL;
  };
  HAWKIimage(string dir,int i) : DbImage(dir), idx(i) {
    objectmaskpt = NULL;
    skymappt = NULL;
    nimusept = NULL;
    skysubtpt = NULL;
  };
  ~HAWKIimage() {
    if (objectmaskpt != NULL) delete(objectmaskpt);
    if (skymappt != NULL) delete(skymappt);
    if (nimusept != NULL) delete(nimusept);
    if (skysubtpt != NULL) delete(skysubtpt);
  };

  int flatFieldImage();

  const string skymapname () { return string("skymap"); }
  const string skyrmsname () { return string("skyrms"); }
  const string nimusename () { return string("nimuse"); }

  FITS* objectmask() {
    getImage(&objectmaskpt,"objectmask");
    return objectmaskpt;
  }
  FITS* objectmask(int ext) {
    getImage(&objectmaskpt,"objectmask",ext);
    return objectmaskpt;
  }

  FITS* skymap() {
    getImage(&skymappt,"skymap");
    return skymappt;
  }
  FITS* skymap(int ext) {
    getImage(&skymappt,"skymap",ext);
    return skymappt;
  }

  FITS* skysubtracted() {
    getImage(&skysubtpt,"skysubtracted");
    return skysubtpt;
  }
  FITS* skysubtracted(int ext) {
    getImage(&skysubtpt,"skysubtracted",ext);
    return skysubtpt;
  }

  // Numer of images used to build the sky value for each
  // pixel.
  //
  FITS* nimuse() {
    getImage(&nimusept,"nimuse");
    return nimusept;
  }
  FITS* nimuse(int ext) {
    getImage(&nimusept,"nimuse",ext);
    return nimusept;
  }

};

#endif
