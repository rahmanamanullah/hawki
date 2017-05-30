#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
#include <list>
#include <sys/stat.h>

// Handle input arguments
#include <unistd.h> 
#include <getopt.h>
#include <string.h>

#include "hawki_image.h"
#include "stats.h"

using namespace std;

bool runccode = true;

const string fitsext(".fits");


// Option and argument parsing
//
static const char *optString = "r:u:m:h";

static const struct option longOpts[] = {
    { "nrej", required_argument, NULL, 'r' },
    { "nuse", required_argument, NULL, 'u' },
    { "nmin", required_argument, NULL, 'm' },
    { "help", no_argument, NULL, 'h' },
    { NULL, no_argument, NULL, 0 }
};


const int NHDU = 4;

typedef HAWKIimage* HAWKIimage_pt;
typedef vector<HAWKIimage_pt> HAWKIlist;
typedef vector<HAWKIimage_pt>::iterator HAWKIlistIterator;

vector<string> readdblist(char *filename);
int imidxinvector (string dbimage, vector<string> dbimages);
bool FileExists(string strFilename);



void display_usage( void )
{
    puts( "hawki_buildskymap - build NIR skymap of a series of images" );
    puts( "" );
    puts( "\thawki_buildskymap [--nres|--nuse|--nmin] <dbtextfile> <dbimage> ..." );

    /* ... */
    exit( EXIT_FAILURE );
}

// - First argument is a text file with all DbImages in chronological
//   order that will be used for building the sky subtraction.
// - Following arguments are DbImages (that also are listed in the text
//   file) for which the skyframe will be constructed.
//
int main(int argc, char *argv[]) {

  // Initialize
  //
  int nrej = 3, nuse = 7, nmin = 10; // Parameters used by ESO pipeline
  int opt = 0;
  int longIndex = 0;

  // Handle input arguments
  //
  //  opterr = 0;

  opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
  while( opt != -1 ) {
    switch( opt ) {
    case 'r' :
      nrej = atoi(optarg);
      break;
    case 'u' :
      nuse = atoi(optarg);
      break;
    case 'm' :
      nmin = atoi(optarg);
      break;
    case 'h':   /* fall-through is intentional */
    case '?':
      display_usage();
    break;
    
    default:
      /* You won't actually get here. */
      break;
    }
    
    opt = getopt( argc, argv, optString );

    //    cout << opt << endl;
  }
    
  char **inputFiles = argv + optind;
  int numInputFiles = argc - optind;

  char *dbfilename = inputFiles[0];
  vector<string> dbimages = readdblist(dbfilename);
  float skyval, skymean, skystd;

  const bool forcefirstpass = false;
  
  // Remove any existing files.
  // 
  for (int dbinr=1; dbinr<numInputFiles; dbinr++) {
    string dir = string(inputFiles[dbinr]);
    HAWKIimage dbimage(dir);

    // Remove old output images if they already
    // exists
    //
    string skymapfile = dbimage.getdbimdir() +
      string("/") + dbimage.skymapname() + fitsext;
    remove(skymapfile.c_str());  // Delete old file if it already exists

    string skyrmsfile = dbimage.getdbimdir() +
      string("/") + dbimage.skyrmsname() + fitsext;
    remove(skyrmsfile.c_str());
    
    string nuseimfile = dbimage.getdbimdir() +
      string("/") + dbimage.nimusename() + fitsext;
    remove(nuseimfile.c_str());
  }

  // Build sky image for one chip at a time to save memory!
  //
  for (int ext=1; ext<=NHDU; ext++) {
    cout << "FITS HDU: " << ext << endl << endl;
    
    int v1 = NULL, v2 = NULL, firstIdx = NULL;
    typedef valarray<float> floatvector;
    typedef valarray<unsigned int> usintgvector;
    vector<floatvector> imagpxBefore, imagpxAfter;
    vector<floatvector> maskpxBefore, maskpxAfter;
    valarray<float> currentImag, currentMask;
    vector<float> mediansky;
    string extname;

    // Go through all DbImages
    //
    for (int dbinr=1; dbinr<numInputFiles; dbinr++) {

      // Name of the DbImage
      //
      string dbimg = string(inputFiles[dbinr]);
      dbimg = dbimg.substr(0,dbimg.find_last_not_of("/")+1);  // Trailing /

      // Determine range of images to use for building the
      // current skymap
      //
      int idx = imidxinvector(dbimg,dbimages);  // current frame
      int i1 = idx-nuse;                        // start frame
      int i2 = idx+nuse;                        // end frame
      if (i1 < 0 ) i1 = 0;
      if (i2 >= int(dbimages.size())) i2 = dbimages.size() - 1;

      if (dbinr == 1) {          // Initialize indices for first iteration
	firstIdx = i1;
	v1 = i1;
	v2 = i1;
      }
 
      // UPDATE FRAME STACKS
      //
      // We work with two stacks, one for the images obtained before
      // the current image and one for the images obtained after.  To limit
      // disk IO, we want to keep frames that have been loaded and will be
      // used again.  First we erase frames from the before stack, and then
      // add images to the after stack.
      //
      while ( v1 != v2 && v1 < i1 && imagpxBefore.size() > 0) {
	imagpxBefore.erase(imagpxBefore.begin());
	maskpxBefore.erase(maskpxBefore.begin());
	v1++;
      }
      while ( v1 == v2 || v2 <= i2 ) {
	HAWKIimage hawkiimage(dbimages[v2]);
	
	FITS* f = hawkiimage.flatfielded(ext);
	FITS* b = hawkiimage.bpm(ext);

	ExtHDU& fHDU = f->extension(ext);
	ExtHDU& bHDU = b->extension(ext);
	
	// We need the extension name later
	//
	if (extname.size() == 0) extname.assign(fHDU.name());

	cout << "\tLOAD: " << dbimages[v2] << "[" << ext << "]" << endl;

	valarray<float> image, bpm;
	valarray<unsigned int> obj;
	fHDU.read(image);
	bHDU.read(bpm);
	valarray<float> mask(bpm);

	/*
	cout << "\tMASK (BPM): " << 100.0*(mask.sum()/mask.size()) <<
	  " % masked" << endl;
	*/

	// Load object mask if it exists
	//
	string objectmaskfile = hawkiimage.getdbimdir() +
	  string("/objectmask.fits");
	if (!forcefirstpass && FileExists(objectmaskfile)) {
	  FITS* o = hawkiimage.objectmask();
	  ExtHDU& oHDU = o->extension(ext);
	  valarray<float> omask;
	  oHDU.read(omask);

	  for (unsigned int i=0; i<omask.size(); i++) {
	    if (omask[i]) mask[i] = 1;              // Add to mask
	  }

	  /*
	  cout << "\tMASK (OBJECTS): " << 100.0*(omask.sum()/omask.size()) <<
	    " % masked" << endl;
	  */
	}

	// Clipped statatistics to find out the sky value.  We are using
	// quite agressive clipping here which will bias the RMS.
	//
	clippedStats<3,4>(image,mask,skyval,skymean,skystd);

	// The object mask will not mask cosmics or moving objects
	// (asteroids, satellits), so we also apply a sigma cut
	//
	// Not sure this is a good idea, the sky varies over the chip
	// so skyval could be off the actual sky value for a given pixel
	// and, as was pointed out above, the RMS is biased.
	//
	unsigned int maskskysigma = 0;
	if (maskskysigma>0) {
	  unsigned int nskymasked = 0;
	  for (unsigned int i=0; i<image.size(); i++) {
	    if ( mask[i] == 0 && fabs(image[i]-skyval) > maskskysigma*skystd ) {
	      mask[i] = 1;
	      nskymasked++;
	    }
	  }

	  cout << "\tMASK (BRIGHT PIX): " <<
	    100.0*((float) nskymasked/(float) image.size()) <<
	    " % masked" << endl;
	}

	cout << "\tMASK (TOTAL): " << 100.0*(mask.sum()/mask.size()) <<
	  " % masked" << endl;
	cout << "\tSKY VAL = " << skyval << " RMS = " << skystd << endl;

	// Normalize with the median sky.  The idea here is to assume that
	// the spatial sky variations are much slower than the temporal.  In
	// other words we use the 'nuse' frames before and after to map the
	// spatial variations and then normalize the map to the the median
	// sky value of the current frame.
	//
	image /= skyval;
	mediansky.push_back((float) skyval);   // We need this later for renormalizing

	// Add image and mask to stack, if we have loaded the image
	// for which we wish to build the stack, save it as "current".
	//
	if (v2 < idx) {
	  imagpxBefore.push_back(image);
	  maskpxBefore.push_back(mask);
	} else if (v2 == idx) {
	  currentImag.resize(image.size());
	  currentImag = image;
	  currentMask.resize(mask.size());
	  currentMask = mask;
	} else {
	  imagpxAfter.push_back(image);
	  maskpxAfter.push_back(mask);
	}
	  
	v2++;
      }

      cout << "\tCURRENT : " << dbimg << " " << idx <<
	" : " << i1 << " - " << i2 << " : " <<
	"sky = " << mediansky[idx-firstIdx] << endl;
     
      const int npix = currentImag.size();      // Total number of pixels

      // Number or images used for creating sky map of
      // current image
      //
      const int nim = imagpxBefore.size() + imagpxAfter.size();

      cout << "\tUsing " << nim << " frames to build sky" << endl;
      /*
      for (int i=0; i<imagpxBefore.size(); i++)
	cout << "\t" << imagpxBefore[i] << endl;
      for (int i=0; i<imagpxAfter.size(); i++)
	cout << "\t" << imagpxAfter[i] << endl;
      */

      valarray<float> skymap(npix), skyrms(npix), iratio(npix);
      valarray<unsigned int> nimuse(npix);
      for (int n=0; n<npix; n++) {

	// Image stack column for given pixel
	//
	valarray<float> skycol(nim);
	valarray<float> mskcol(nim);
	for (unsigned int i=0; i<imagpxBefore.size(); i++) {
	  skycol[i] = imagpxBefore[i][n];
	  mskcol[i] = maskpxBefore[i][n];
	}
	for (unsigned int i=0; i<imagpxAfter.size(); i++) {
	  skycol[imagpxBefore.size()+i] = imagpxAfter[i][n];
	  mskcol[maskpxBefore.size()+i] = maskpxAfter[i][n];
	}

	// Calculate sky value for given pixel, if there are not at
	// least nmin available, we use the median.
	//
	float median, mean, rms, nused;
	
	runningStats(skycol,mskcol,nrej,median,mean,rms,nused);
	float skyinpix = mean;
	if (nused < nmin) {  // Fall back on no rejection
	  runningStats<0>(skycol,mskcol,median,mean,rms,nused);
	  skyinpix = median;
	}

	skymap[n] = nused > 0 ? skyinpix : 0.0;
	skyrms[n] = nused > 0 ? rms : -1.0;
	nimuse[n] = int(nused);

	//	if (idx == 50)
	//	  cout << mediansky.size() << " " << idx-firstIdx << " " << i1 << " " << i2 << " " << currentImag[n] << " " << mediansky[idx-firstIdx] << endl;
		
	// Calculate the ratio between the image for which the
	// skymap is being build and the sky map.  This will be
	// used later to rescale the the map.
	//
	if (skymap[n] > 0.0)
	  iratio[n] = mediansky[idx-firstIdx]*currentImag[n]/skymap[n];
	else
	  iratio[n] = 0.0;
      }

      // Calculate median ratio between image and sky
      //
      clippedStats<3,4>(iratio,currentMask,skyval,skymean,skystd);
      cout << "\tmedian image/sky ratio = " << skyval << endl;
      skymap *= skyval;
      skyrms *= skyval;

      // Save the skymap, skyrms and nimuse to FITS files
      //
      HAWKIimage dbimage(dbimg);

      string skymapfile =  dbimage.getdbimdir() + string("/") +
	dbimage.skymapname() + fitsext;
      FITS skymapfits(skymapfile,Write,false);
      if (ext == 1) {
	PHDU& phdu = skymapfits.pHDU();
	phdu.addKey("NREJ",nrej,"Number of high and low rejected");
	phdu.addKey("NUSE",nuse,"Number of frames before and after used");
	phdu.addKey("NMIN",nmin,"Min number of frames for building sky");
      }
      
      string skyrmsfile =  dbimage.getdbimdir() + string("/") +
	dbimage.skyrmsname() + fitsext;
      FITS skyrmsfits(skyrmsfile,Write,false);

      string nimusefile =  dbimage.getdbimdir() + string("/") +
	dbimage.nimusename() + fitsext;
      FITS nimusefits(nimusefile,Write,false);
      
      vector<long> extAx(2,2048);
      long fpixel(1);

      ExtHDU* skymapext = skymapfits.addImage(extname,FLOAT_IMG,extAx);
      skymapext->addKey("ORIGSKY",mediansky[idx-firstIdx],"Median sky");
      skymapext->addKey("SKYRATI",skyval,"Ratio between image and sky map");
      skymapext->write(fpixel,npix,skymap);

      ExtHDU* skyrmsext = skyrmsfits.addImage(extname,FLOAT_IMG,extAx);
      skyrmsext->write(fpixel,npix,skyrms);

      ExtHDU* nimuseext = nimusefits.addImage(extname,SHORT_IMG,extAx);
      nimuseext->write(fpixel,npix,nimuse);
      
      // Update image stacks
      //
      imagpxBefore.push_back(currentImag);
      maskpxBefore.push_back(currentMask);

      if (imagpxAfter.size() > 0) {
	currentImag.resize(imagpxAfter[0].size());
	currentImag = imagpxAfter[0];
	currentMask.resize(maskpxAfter[0].size());
	currentMask = maskpxAfter[0];
      
	imagpxAfter.erase(imagpxAfter.begin());
	maskpxAfter.erase(maskpxAfter.begin());
      }

      cout << endl;
    }
  }
}

// Find index of dbimage in vector
//
int imidxinvector (string dbimage, vector<string> dbimages) {
  for (unsigned int i=0; i<dbimages.size(); i++) {
    string::size_type loc = dbimages[i].find(dbimage,0);
    if ( loc != string::npos)
      return i;
  }

  return -1;
}

// Read the DbImage list from text file
//
vector<string> readdblist(char *filename) {
  ifstream dbframes(filename);

  const int max_num_of_char_in_a_line = 256;
  vector<string> dbimages;
  
  while(!dbframes.eof()) {
    char dbimage[max_num_of_char_in_a_line+1];
    dbframes.getline(dbimage,max_num_of_char_in_a_line);
    if ( dbimage[0] != '\n' && dbimage[0] != '#' && dbimage[0] != '\0' )
      dbimages.push_back(dbimage);
  }

  dbframes.close();

  return dbimages;
}

bool FileExists(string strFilename) {
  struct stat stFileInfo;
  bool blnReturn;
  int intStat;

  // Attempt to get the file attributes
  intStat = stat(strFilename.c_str(),&stFileInfo);
  if(intStat == 0) {
    blnReturn = true;
  } else {
    blnReturn = false;
  }
  
  return(blnReturn);
}

// EOF
//

