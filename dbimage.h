#ifndef DBIMAGE_H
#define DBIMAGE_H

#include <sys/stat.h> 
#include <string>
#include <iostream>
#include <CCfits/CCfits>
#include <cmath>

using namespace std;
using namespace CCfits;

class DbImage {
 private :
  string dbimdir;
  FITS *rawpt, *darkpt, *bpmpt, *deadpt, *flatpt, *flatfieldedpt;

 protected :
  bool fileexists(string strFilename);

  FITS* readFile(string type);               // Read the full file
  FITS* readFile(string type,int ext);       // Read only a single extension

  void getImage(FITS **fitspt,string type) {
    if (*fitspt == NULL) *fitspt = readFile(type);    // Read if needed
  }
  void getImage(FITS **fitspt,string type, int ext) { // Get a given extension
    if (*fitspt == NULL) *fitspt = readFile(type,ext);
  }

  int findfile(string* type);

  const string flatfieldedname () { return string("flatfielded"); }

 public :
  DbImage() { DbImage(string("0")); };
  DbImage(string dir) : dbimdir(dir), rawpt(NULL), darkpt(NULL),
    bpmpt(NULL), deadpt(NULL), flatpt(NULL), flatfieldedpt(NULL)
    {
      FITS::setVerboseMode(false);
    }
  ~DbImage() {
    if (rawpt  != NULL) delete(rawpt);
    if (bpmpt  != NULL) delete(bpmpt);
    if (deadpt != NULL) delete(deadpt);
    if (darkpt != NULL) delete(darkpt);
    if (flatpt != NULL) delete(flatpt);
    if (flatfieldedpt != NULL) delete(flatfieldedpt);
  }

  /*
  DbImage (const DbImage &dbi) {
    dbimdir = dbi.dbimdir;
    rawpt  = dbi.rawpt;
    darkpt = dbi.darkpt;
    bpmpt  = dbi.bpmpt;
    deadpt = dbi.deadpt;
    flatpt = dbi.flatpt;
  }
  DbImage operator=(const DbImage &dbi) {
    dbimdir = dbi.dbimdir;
    rawpt  = dbi.rawpt;
    darkpt = dbi.darkpt;
    bpmpt  = dbi.bpmpt;
    deadpt = dbi.deadpt;
    flatpt = dbi.flatpt;

    return *this;
  }
  */

  void setdbimdir(string dir) {this->dbimdir = dir;}
  string getdbimdir() {return this->dbimdir + string("/");}
 
  FITS* dark() {
    getImage(&darkpt,"dark");
    return darkpt;
  }
  FITS* dark(int ext) {
    getImage(&darkpt,"dark",ext);
    return darkpt;
  }

  FITS* flat() {
    getImage(&flatpt,"flatcube");
    return flatpt;
  }
  FITS* flat(int ext) {
    getImage(&flatpt,"flatcube",ext);
    return flatpt;
  }

  FITS* raw() {
    getImage(&rawpt,"raw");
    return rawpt;
  }
  FITS* raw(int ext) {
    getImage(&rawpt,"raw",ext);
    return rawpt;
  }

  FITS* bpm() {
    getImage(&bpmpt,"bpm");
    return bpmpt;
  }
  FITS* bpm(int ext) {
    getImage(&bpmpt,"bpm",ext);
    return bpmpt;
  }

  FITS* dead() {
    getImage(&deadpt,"dead");
    return deadpt;
  }
  FITS* dead(int ext) {
    getImage(&deadpt,"dead",ext);
    return deadpt;
  }

  FITS* flatfielded() {
    getImage(&flatfieldedpt,flatfieldedname());
    return flatfieldedpt;
  }
  FITS* flatfielded(int ext) {
    getImage(&flatfieldedpt,flatfieldedname(),ext);
    return flatfieldedpt;
  }


  // virtual int flatfieldimage();
};

typedef auto_ptr<DbImage> DbImagePt;
typedef vector<DbImagePt> DbImageList;

#endif
