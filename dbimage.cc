#include "dbimage.h"

FITS* DbImage::readFile(string type) {
  string file = type;

  FITS* pInfile = NULL;
  if (findfile(&file)) {
    cerr << "Warning: no " << type << " file!" << endl;
  } else {
    try {
      pInfile = new FITS(file,Read,true);
    } catch (FITS::CantOpen) {
      cout << "ERROR: failed to open: "+file << endl;
    }
  }

  return pInfile;
}

FITS* DbImage::readFile(string type,int ext) {
  string file = type;

  FITS* pInfile = NULL;
  if (findfile(&file)) {
    cerr << "Warning: no " << type << " file!" << endl;
  } else {
    try {
      pInfile = new FITS(file,Read,ext,true);
    } catch (FITS::CantOpen) {
      cout << "ERROR: failed to open: "+file << endl;
    }
  }

  return pInfile;
}

///////////////////////////////////////////////////////////////////////
//
// PROTECTED


bool DbImage::fileexists(string strFilename) {
  struct stat stFileInfo;
  bool blnReturn;
  int intStat;

  // Attempt to get the file attributes
  intStat = stat(strFilename.c_str(),&stFileInfo);
  if(intStat == 0) {
    // We were able to get the file attributes
    // so the file obviously exists.
    blnReturn = true;
  } else {
    // We were not able to get the file attributes.
    // This may mean that we don't have permission to
    // access the folder which contains this file. If you
    // need to do that level of checking, lookup the
    // return values of stat which will give you
    // more details on why stat failed.
    blnReturn = false;
  }
  
  return(blnReturn);
}


int DbImage::findfile(string* type) {
  string file = getdbimdir() + "/" + *type + ".fits";

  if ( ! fileexists(file) )
    file = file + ".gz";
  if ( fileexists(file) ) {
    *type = file;
    return 0;
  } else {
    return -1;
  }
}
